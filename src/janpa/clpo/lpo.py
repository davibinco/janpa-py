# janpa/clpo/lpo.py
"""LPO / CLPO (Localized / Chemist's Localized Property-optimized) orbitals. Replaces CLPO/PropertyOptimizedOrbitals.java. """
from __future__ import annotations
from dataclasses import dataclass, field
import re
import numpy as np
import networkx as nx

from janpa.gto.basis import AtomicCenter, BasisFunction
from janpa.utils.matrix import (
    transform_symmetric_to_new_basis,
    symmetric_orthogonalize,
)
from janpa.utils.simult_diag import simultaneous_diagonalization
from .lodesc import LODescription, LO_TYPE_RY, LO_TYPE_LP, LO_TYPE_BD, LO_TYPE_NB
from .hybrids import AtomicHybrids

@dataclass
class CLPOOptions:
    """Options for LPO/CLPO computation."""
    hybr_opt_conv_thresh: float = 1e-5
    hybr_opt_max_iter: int = 10000
    max_bond_ionicity_threshold: float = 0.90
    ry_occ_print_threshold: float = 1e-3
    edges: str = ""

@dataclass
class CLPOResult:
    """Container for LPO/CLPO results."""
    lpo: LODescription | None = None
    clpo: LODescription | None = None
    sds_in_hybrid_basis: np.ndarray | None = None

class PropertyOptimizedOrbitals:
    def __init__(self):
        self.sds_nao = None
        self.n_atoms = 0
        self.naos_at_center = None
        self.n_naos = 0
        self.hybrids_of_atoms = None
        self.dab = None
        self.lpo_descript = None
        self.clpo_descript = None
        self.options = None
        self.centers = None
        self.target_function_conv_thresh = 1e-5
        self.opt_max_iter = 10000
        self.opt_lewis_mode = False
        self.optimize_hybrids_n_iters_done = -1

    def get_naos_at_centers(self, naos):
        counts = np.zeros(self.n_atoms, dtype=int)
        for b in range(len(naos)):
            counts[naos[b].center_id - 1] += 1
        
        result = [np.zeros(c, dtype=int) for c in counts]
        counters = np.zeros(self.n_atoms, dtype=int)
        
        for b in range(len(naos)):
            c = naos[b].center_id - 1
            result[c][counters[c]] = b
            counters[c] += 1
            
        return result

    def diatomic_submatrices(self):
        result = []
        for i in range(self.n_atoms):
            row = []
            for j in range(i + 1):
                row.append(self.sds_nao[np.ix_(self.naos_at_center[i], self.naos_at_center[j])])
            result.append(row)
        return result

    def dihybrid_matrix(self, atom_a, i_a, atom_b, i_b):
        daah = self.dab_dot_h(atom_a, atom_a, i_a)
        dbbh = self.dab_dot_h(atom_b, atom_b, i_b)
        dab_h = self.dab_dot_h(atom_a, atom_b, i_b)
        d11 = self.hybrids_of_atoms[atom_a].hybrid_scalar_mul(i_a, daah)
        d22 = self.hybrids_of_atoms[atom_b].hybrid_scalar_mul(i_b, dbbh)
        d12 = self.hybrids_of_atoms[atom_a].hybrid_scalar_mul(i_a, dab_h)
        result = np.zeros((2, 2))
        result[0, 0] = d11
        result[0, 1] = d12
        result[1, 0] = d12
        result[1, 1] = d22
        return result

    def _print_hybrid_stat_table(self, skip_pairs, skip_singles, occ_threshold):
        print(" AtomA (VecId) AtomB (VecId) occupancy Dij")
        occ_sum_above_thresh = 0.0
        occ_sum_below_thresh = 0.0
        n_items_above_thresh = 0
        n_items_below_thresh = 0
        for a in range(self.n_atoms):
            for h in range(self.hybrids_of_atoms[a].n_valid_hybrids):
                h_a = self.hybrids_of_atoms[a].get_hybrid(h)
                b = self.hybrids_of_atoms[a].friend_atom_index[h]
                i_b = self.hybrids_of_atoms[a].friend_hybrid_index[h]
                if skip_pairs and (b != -1): continue
                if skip_singles and (b == -1): continue
                dii = float((h_a.T @ self.dab[a][a] @ h_a)[0, 0])
                if dii < occ_threshold:
                    occ_sum_below_thresh += dii
                    n_items_below_thresh += 1
                    continue
                else:
                    n_items_above_thresh += 1
                    occ_sum_above_thresh += dii
                
                b_str = f"{b+1:5d}" if b != -1 else f"{'-1':>5s}"
                ib_str = f"{i_b+1:5d}" if i_b != -1 else f"{'-1':>5s}"
                print(f" {a+1:5d} ({h+1:5d}) {b_str} ({ib_str}) {dii:9.5f} ", end="")
                if b != -1:
                    val = float((h_a.T @ self.dab_dot_h(a, b, i_b))[0, 0])
                    print(f"{val:7.3f}")
                else:
                    print(f"{'n/a':>7s}")
                    
        if occ_threshold > 0:
            print(f"Total: {occ_sum_above_thresh + occ_sum_below_thresh:.3f} electrons = {occ_sum_above_thresh:.3f} electrons in {n_items_above_thresh} printed items ")
            print(f"and {occ_sum_below_thresh:.3f} electrons in {n_items_below_thresh} items (each having occupancy below threshold = {occ_threshold:.4f})")
        else:
            print(f"Total: {occ_sum_above_thresh:.3f} electrons in {n_items_above_thresh} printed items ")
        print()

    def dab_dot_h(self, a, b, i_b):
        if b <= a:
            _dij = self.dab[a][b]
            u_col = self.hybrids_of_atoms[b].u[:, i_b]
            result = _dij @ u_col
        else:
            _dij = self.dab[b][a]
            u_col = self.hybrids_of_atoms[b].u[:, i_b]
            result = _dij.T @ u_col
        return result.reshape(-1, 1)

    def grad_dnlo_sq(self, a, pwin):
        u_a = self.hybrids_of_atoms[a].u
        n_vecs = u_a.shape[1]
        vec_len = u_a.shape[0]
        result = np.zeros((vec_len, n_vecs))
        win = 0.0
        for i in range(n_vecs):
            daa_ha = self.dab_dot_h(a, a, i)
            wii = self.hybrids_of_atoms[a].hybrid_scalar_mul(i, daa_ha)
            b = self.hybrids_of_atoms[a].friend_atom_index[i]
            if b == -1:
                result[:, i] = daa_ha[:, 0] * wii
                win += wii * wii
            else:
                dab_hb = self.dab_dot_h(a, b, self.hybrids_of_atoms[a].friend_hybrid_index[i])
                wij = self.hybrids_of_atoms[a].hybrid_scalar_mul(i, dab_hb)
                result[:, i] = daa_ha[:, 0] * wii + dab_hb[:, 0] * wij
                win += wii * wii + wij * wij
        if pwin is not None:
            pwin[0] = win
        return result

    def grad_dnlo_occ2(self, a, pwin):
        u_a = self.hybrids_of_atoms[a].u
        n_vecs = u_a.shape[1]
        vec_len = u_a.shape[0]
        result = np.zeros((vec_len, n_vecs))
        win = 0.0
        for i in range(n_vecs):
            daa_ha = self.dab_dot_h(a, a, i)
            dii = self.hybrids_of_atoms[a].hybrid_scalar_mul(i, daa_ha)
            b = self.hybrids_of_atoms[a].friend_atom_index[i]
            if b == -1:
                wii = dii
                result[:, i] += daa_ha[:, 0] * wii
                win += dii * dii
            else:
                hb = self.hybrids_of_atoms[a].friend_hybrid_index[i]
                dab_hb = self.dab_dot_h(a, b, hb)
                dij = self.hybrids_of_atoms[a].hybrid_scalar_mul(i, dab_hb)
                dbb_hb = self.dab_dot_h(b, b, hb)
                djj = self.hybrids_of_atoms[b].hybrid_scalar_mul(hb, dbb_hb)
                
                sqrt_val = np.sqrt((dii - djj)**2 + 4 * dij**2)
                bd_occ = 0.5 * (dii + djj + sqrt_val)
                wii = bd_occ * 0.5 * (1 + (dii - djj) / sqrt_val)
                wij = bd_occ * dij / sqrt_val
                
                result[:, i] += daa_ha[:, 0] * wii + dab_hb[:, 0] * wij
                win += bd_occ * bd_occ / 2.0
        if pwin is not None:
            pwin[0] = win
        return result

    def grad_dnlo(self, a, pwin):
        if not self.opt_lewis_mode:
            return self.grad_dnlo_sq(a, pwin)
        else:
            return self.grad_dnlo_occ2(a, pwin)

    @staticmethod
    def parse_edges(input_str):
        edges = []
        start_index = input_str.find('[')
        end_index = input_str.find(']')
        if start_index == -1 or end_index == -1 or start_index >= end_index:
            raise ValueError("Wasn't found a valid list delimited by the brakets '[...]' on the input.")
        content = input_str[start_index + 1:end_index]
        pattern = re.compile(r"\(([^)]*)\)")
        for match in pattern.finditer(content):
            cedge = match.group(1).strip()
            if not cedge:
                raise ValueError("Empty edge found '()', not allowed.")
            if cedge.endswith(","):
                cedge = cedge[:-1].strip()
            parts = cedge.split(",")
            try:
                if len(parts) == 1:
                    node = int(parts[0].strip())
                    edges.append([node, node])
                elif len(parts) == 2:
                    origin = int(parts[0].strip())
                    destiny = int(parts[1].strip())
                    edges.append([origin, destiny])
                else:
                    raise ValueError(f"Edge with invalid length ({len(parts)} elements). Only 1 or 2 allowed. Found: ({match.group(1)})")
            except ValueError:
                raise ValueError(f"Not numerical value found on the edge: ({match.group(1)})")
        return edges

    def contains_edge(self, lst, a, b):
        for e in lst:
            if (e[0] == a and e[1] == b) or (e[0] == b and e[1] == a):
                return True
        return False

    def optimize_hybrids(self):
        iter_count = 0
        g = [None] * self.n_atoms
        for a in range(self.n_atoms):
            self.hybrids_of_atoms[a].backup_u()
            
        daa_h_basis = [None] * self.n_atoms
        for a in range(self.n_atoms):
            daa_h_basis[a] = self.hybrids_of_atoms[a].u.T @ self.dab[a][a] @ self.hybrids_of_atoms[a].u
            
        changes = np.zeros((2, 2))
        pwin = [0.0]
        win_prev = -1.0
        print("Iterative optimization of atomic hybrids")
        print(" Target Total max.atomic Total max.atomic Nxt.step")
        print(" Itr. function new-old ||Unew-U||^2 ||unew-u||^2 ||Dnew-D||^2 ||dnew-d||^2 lambda ")
        
        invlambda = 0.0
        dont_scale_stepsize = True
        
        while True:
            iter_count += 1
            self.optimize_hybrids_n_iters_done = iter_count
            win = 0.0
            invlambda0 = 0.0
            dfdlam = 0.0
            
            for a in range(self.n_atoms):
                g[a] = self.grad_dnlo(a, pwin)
                win += pwin[0]
                invlambda0 += np.abs(np.trace(g[a].T @ self.hybrids_of_atoms[a].u))
                dfdlam += np.trace(g[a] @ g[a].T)
                tmp = g[a].T @ self.hybrids_of_atoms[a].u
                dfdlam -= np.trace(tmp @ tmp)
                
            invlambda0 = invlambda0 / self.n_naos
            print(f"{iter_count:5d} {win:12.4f} {(win) if iter_count == 1 else (win - win_prev):+8.1e}", end="")
            
            if win > win_prev:
                invlambda = 1e-16 * invlambda0
                dont_scale_stepsize = True
                for a in range(self.n_atoms):
                    self.hybrids_of_atoms[a].backup_u()
                if win - win_prev < self.target_function_conv_thresh:
                    print(f" << Converged! (threshold = {self.target_function_conv_thresh:.2e}) >>")
                    print("Optimization of hybrids finished")
                    return win
                win_prev = win
            else:
                for a in range(self.n_atoms):
                    self.hybrids_of_atoms[a].restore_u()
                if win_prev - win < self.target_function_conv_thresh:
                    print(" << No step can be taken >>")
                    print("Optimization of hybrids finished")
                    return win_prev
                if dont_scale_stepsize:
                    invlambda = invlambda0
                else:
                    invlambda *= 2.0
                dont_scale_stepsize = False
                
            if iter_count >= self.opt_max_iter:
                print(" << Maximum number of iterations reached >>")
                print("Optimization of hybrids finished")
                for a in range(self.n_atoms):
                    self.hybrids_of_atoms[a].restore_u()
                return win_prev
                
            changes[:] = 0.0
            for a in range(self.n_atoms):
                ua_new = symmetric_orthogonalize(g[a] + self.hybrids_of_atoms[a].u * invlambda)
                du2 = np.linalg.norm(ua_new - self.hybrids_of_atoms[a].u, 'fro')**2
                changes[0, 0] += du2
                n_dim = len(self.hybrids_of_atoms[a].nao_indices)
                if du2 / n_dim > changes[0, 1]:
                    changes[0, 1] = du2 / n_dim
                    
                daa_new = ua_new.T @ self.dab[a][a] @ ua_new
                du2 = np.linalg.norm(daa_h_basis[a] - daa_new, 'fro')**2
                daa_h_basis[a] = daa_new
                changes[1, 0] += du2
                if du2 > changes[1, 1]:
                    changes[1, 1] = du2
                self.hybrids_of_atoms[a].u = ua_new
                
            print(f"{changes[0, 0]:12.6f} {changes[0, 1]:13.7f} {changes[1, 0]:12.5f} {changes[1, 1]:13.7f} ", end="")
            if dont_scale_stepsize:
                print(f"{'(none)':>10s}")
            else:
                print(f"{1/invlambda:10.2e}")

    def get_global_hybrid_indices(self):
        addr2i_h = [None] * self.n_atoms
        i_h = 0
        for a in range(self.n_atoms):
            a_naos = self.hybrids_of_atoms[a].nao_indices
            addr2i_h[a] = np.zeros(len(a_naos), dtype=int)
            for h in range(self.hybrids_of_atoms[a].n_valid_hybrids):
                addr2i_h[a][h] = i_h
                i_h += 1
        return addr2i_h

    def sds_in_hybrid_basis(self, hybr_addresses, nao2ho_holder):
        n_tot_hybrids = sum(self.hybrids_of_atoms[a].n_valid_hybrids for a in range(self.n_atoms))
        d_oho = np.zeros((n_tot_hybrids, n_tot_hybrids))
        
        for a in range(self.n_atoms):
            a_naos = self.hybrids_of_atoms[a].nao_indices
            if nao2ho_holder is not None:
                for h in range(self.hybrids_of_atoms[a].n_valid_hybrids):
                    vec = self.hybrids_of_atoms[a].u[:, h]
                    i_h = hybr_addresses[a][h]
                    nao2ho_holder.nao_to_hybrids[a_naos, i_h] = vec
                    
            daa = self.hybrids_of_atoms[a].u.T @ self.dab[a][a] @ self.hybrids_of_atoms[a].u
            d_oho[np.ix_(a_naos, a_naos)] = daa * 0.5
            
            for b in range(a):
                b_naos = self.hybrids_of_atoms[b].nao_indices
                dab = self.hybrids_of_atoms[a].u.T @ self.dab[a][b] @ self.hybrids_of_atoms[b].u
                d_oho[np.ix_(a_naos, b_naos)] = dab
                
        d_oho = d_oho + d_oho.T
        return d_oho

    def print_bond_matrix(self, num_bonds):
        print("Number of two-center(2C) BD orbitals for each pair of atoms")
        print(f"{'Centr. A/B':>10s}", end="")
        for i in range(self.n_atoms):
            print(f"{i+1:10d}", end="")
        print()
        for i in range(self.n_atoms):
            print(f"{i+1:7d} ", end="")
            for j in range(i):
                print(f"{'':>10s}", end="")
            print(f"({num_bonds[i][i]:7d})", end="")
            for j in range(i + 1, self.n_atoms):
                print(f"{num_bonds[i][j]:10d}", end="")
            print()

    def num_bonds(self):
        result = np.zeros((self.n_atoms, self.n_atoms), dtype=int)
        for a in range(self.n_atoms):
            for ha in range(self.hybrids_of_atoms[a].n_valid_hybrids):
                b = self.hybrids_of_atoms[a].friend_atom_index[ha]
                if b != -1:
                    result[a, b] += 1
        for a in range(self.n_atoms):
            result[a, a] = 0
            for b in range(self.n_atoms):
                if b != a:
                    result[a, a] += result[a, b]
        return result

    def cs_guess(self):
        ui = [None] * self.n_atoms
        for a in range(self.n_atoms):
            dgs = [None] * self.n_atoms
            for b in range(self.n_atoms):
                if b <= a:
                    dgs[b] = self.dab[a][b] @ self.dab[a][b].T
                else:
                    dgs[b] = self.dab[b][a].T @ self.dab[b][a]
            
            ui[a], _ = simultaneous_diagonalization(dgs)
            self.hybrids_of_atoms[a].u = ui[a]
            
        for a in range(self.n_atoms):
            self.hybrids_of_atoms[a].n_valid_hybrids = len(self.hybrids_of_atoms[a].nao_indices)
    def reconnect_hybrids(self, d_in_hybrid_basis, hybr_addresses, allowall):
        max_bond_ionicity_threshold = self.options.max_bond_ionicity_threshold
        graph_table = [f"{'ID':<4s} {'Description':<30s} {'Occupancy':<10s} Composition"]
        
        n_edges_max = self.n_naos * (2 * self.n_naos - 1) // 2
        if self.opt_lewis_mode:
            n_edges_max += self.n_naos
            
        custom = False
        custom_edges = []
        if self.options.edges:
            custom_edges = self.parse_edges(self.options.edges)
            custom = True
            if len(custom_edges) > n_edges_max:
                raise ValueError(f"More edges than expected provided, provided: {len(custom_edges)}, expected: {n_edges_max}")
                
        G = nx.Graph()
        for i in range(2 * self.n_naos):
            G.add_node(i)
            
        nao_owners = [None] * self.n_naos
        for a in range(self.n_atoms):
            for ha in range(len(self.hybrids_of_atoms[a].nao_indices)):
                self.hybrids_of_atoms[a].friend_atom_index[ha] = -1
                self.hybrids_of_atoms[a].friend_hybrid_index[ha] = -1
                i_a = hybr_addresses[a][ha]
                nao_owners[i_a] = [a, ha]
                d_aa = d_in_hybrid_basis[i_a, i_a]
                
                nadded = True
                if custom:
                    if self.contains_edge(custom_edges, i_a, i_a):
                        G.add_edge(i_a, self.n_naos + i_a, weight=1000.0)
                        nadded = False
                        
                if nadded and self.opt_lewis_mode:
                    G.add_edge(i_a, self.n_naos + i_a, weight=d_aa * d_aa)
                    
                for b in range(self.n_atoms):
                    if b == a: continue
                    for hb in range(len(self.hybrids_of_atoms[b].nao_indices)):
                        i_b = hybr_addresses[b][hb]
                        if i_b > i_a: continue
                        if custom:
                            if self.contains_edge(custom_edges, i_a, i_b):
                                G.add_edge(i_b, i_a, weight=1000.0)
                                continue
                                
                        d_bb = d_in_hybrid_basis[i_b, i_b]
                        d_ab = d_in_hybrid_basis[i_a, i_b]
                        bd_occ = 0.5 * (d_aa + d_bb + np.sqrt((d_aa - d_bb)**2 + 4 * d_ab**2))
                        nb_occ = d_aa + d_bb - bd_occ
                        f = (bd_occ > 1.0) and (nb_occ < 1.0) and (np.cos(np.arctan2(2 * d_ab, np.abs(d_aa - d_bb))) < max_bond_ionicity_threshold)
                        
                        if allowall or f:
                            weight = bd_occ * bd_occ if self.opt_lewis_mode else d_ab * d_ab
                            G.add_edge(i_a, i_b, weight=weight)
                            
        matching = nx.max_weight_matching(G, maxcardinality=False, weight='weight')
        
        remap = np.full(2 * self.n_naos, -1, dtype=int)
        for u, v in matching:
            remap[u] = v
            remap[v] = u
            
        result = 0.0
        row_id = 1
        for i in range(self.n_naos):
            j = remap[i]
            if j >= self.n_naos or j == -1:
                occ = d_in_hybrid_basis[i, i]
                a, ha = nao_owners[i]
                graph_table.append(f"{row_id:4d} (LP) {self.centers[a].name + str(a + 1):<26s} {occ:10.5f} 1.0 * h{ha + 1}@{self.centers[a].name}{a + 1}")
                result += occ * occ
                row_id += 1
                continue
                
            if i < j:
                d11 = d_in_hybrid_basis[i, i]
                d12 = d_in_hybrid_basis[i, j]
                d22 = d_in_hybrid_basis[j, j]
                a, ha = nao_owners[i]
                b, hb = nao_owners[j]
                
                self.hybrids_of_atoms[a].friend_atom_index[ha] = b
                self.hybrids_of_atoms[a].friend_hybrid_index[ha] = hb
                self.hybrids_of_atoms[b].friend_atom_index[hb] = a
                self.hybrids_of_atoms[b].friend_hybrid_index[hb] = ha
                
                bd_occ = 0.5 * (d11 + d22 + np.sqrt((d11 - d22)**2 + 4 * d12**2))
                nb_occ = d11 + d22 - bd_occ
                io = np.cos(np.arctan2(2 * d12, np.abs(d11 - d22)))
                theta = 0.5 * np.arctan2(2 * d12, d11 - d22)
                c = np.cos(theta)
                s = np.sin(theta)
                
                atom_a = f"{self.centers[a].name}{a + 1}"
                atom_b = f"{self.centers[b].name}{b + 1}"
                
                graph_table.append(f"{row_id:4d} (BD) {atom_a}-{atom_b}, Io = {io:.4f} {bd_occ:10.5f} h{ha + 1}@{atom_a} * ({c:7.4f}) + h{hb + 1}@{atom_b} * ({-s:7.4f})")
                row_id += 1
                graph_table.append(f"{row_id:4d} {atom_a}-{atom_b}, antibonding (NB) {nb_occ:10.5f} h{ha + 1}@{atom_a} * ({-s:7.4f}) + h{hb + 1}@{atom_b} * ({-c:7.4f})")
                row_id += 1
                
                result += bd_occ * bd_occ if self.opt_lewis_mode else d11 * d11 + 2 * d12 * d12 + d22 * d22
            
        return result

    def create_los(self, lo_name_str, lo_info):
        addr2i_h = self.get_global_hybrid_indices()
        d_aho = self.sds_in_hybrid_basis(addr2i_h, lo_info)
        i_lpo = 0
        is_lewis = np.zeros(self.n_naos, dtype=bool)
        ry_print_threshold = self.options.ry_occ_print_threshold
        
        print()
        print(f"*** Summary of {lo_name_str} results")
        print()
        print(f" {lo_name_str:>5s}\t{' D e s c r i p t i o n ':^35s}\t{'Occupancy':>9s}\tComposition")
        
        l2_sum = 0.0
        lp_sum = 0.0
        num_lp = 0
        ry_sum = 0.0
        num_ry = 0
        bd_sum = 0.0
        num_bd = 0
        nb_sum = 0.0
        ry_below_thresh_occs = 0.0
        num_ry_below_thresh = 0
        
        for a in range(self.n_atoms):
            for ha in range(self.hybrids_of_atoms[a].n_valid_hybrids):
                i_ha = addr2i_h[a][ha]
                if lo_info is not None:
                    lo_info.host_atom_of_hybrid[i_ha] = a
                    
                b = self.hybrids_of_atoms[a].friend_atom_index[ha]
                if b == -1:
                    occ = d_aho[i_ha, i_ha]
                    is_lewis[i_lpo] = occ > 1.0
                    if lo_info is not None:
                        lo_info.lo_to_hybrids[i_lpo, i_ha] = 1.0
                        lo_info.hybrids_of_lo[i_lpo] = [i_ha]
                    i_lpo += 1
                    
                    host_info = f"{self.centers[a].name}{a + 1}"
                    composition = f"1.0 * h{i_ha + 1}@{host_info}"
                    
                    if occ > 1.0:
                        lp_sum += occ
                        l2_sum += occ * occ
                        num_lp += 1
                        print(f"{i_lpo:5d}\t{'(LP) ' + host_info:<35s}\t{occ:9.5f}\t{composition}")
                        if lo_info is not None:
                            lo_info.lo_labels[i_lpo - 1] = f"{host_info}:LP"
                            lo_info.lo_types[i_lpo - 1] = LO_TYPE_LP
                    else:
                        ry_sum += occ
                        num_ry += 1
                        if lo_info is not None:
                            lo_info.lo_labels[i_lpo - 1] = f"{host_info}:RY"
                            lo_info.lo_types[i_lpo - 1] = LO_TYPE_RY
                        if occ > ry_print_threshold:
                            print(f"{i_lpo:5d}\t{host_info + ' (RY)':^35s}\t{occ:9.5f}\t{composition}")
                        else:
                            ry_below_thresh_occs += occ
                            num_ry_below_thresh += 1
                else:
                    hb = self.hybrids_of_atoms[a].friend_hybrid_index[ha]
                    i_hb = addr2i_h[b][hb]
                    if i_hb > i_ha:
                        d11 = d_aho[i_ha, i_ha]
                        d12 = d_aho[i_ha, i_hb]
                        d22 = d_aho[i_hb, i_hb]
                        d = np.array([[d11, d12], [d12, d22]])
                        vals, vecs = np.linalg.eigh(d)
                        idx = np.argsort(vals)[::-1]
                        vals = vals[idx]
                        cf = vecs[:, idx]
                        
                        bd_sum += vals[0]
                        nb_sum += vals[1]
                        
                        host_atoms = f"{self.centers[a].name}{a + 1}-{self.centers[b].name}{b + 1}"
                        if lo_info is not None:
                            lo_info.lo_to_hybrids[i_lpo, i_ha] = cf[0, 0]
                            lo_info.lo_to_hybrids[i_lpo, i_hb] = cf[1, 0]
                            lo_info.hybrids_of_lo[i_lpo] = [i_ha, i_hb]
                            lo_info.lo_types[i_lpo] = LO_TYPE_BD
                            lo_info.lo_labels[i_lpo] = f"{host_atoms}:BD"
                            
                        is_lewis[i_lpo] = True
                        i_lpo += 1
                        
                        if lo_info is not None:
                            lo_info.lo_to_hybrids[i_lpo, i_ha] = cf[0, 1]
                            lo_info.lo_to_hybrids[i_lpo, i_hb] = cf[1, 1]
                            lo_info.hybrids_of_lo[i_lpo] = [i_ha, i_hb]
                            lo_info.lo_types[i_lpo] = LO_TYPE_NB
                            lo_info.lo_labels[i_lpo] = f"{host_atoms}:NB"
                        i_lpo += 1
                        
                        info1 = f"h{i_ha + 1}@{self.centers[a].name}{a + 1} * ({cf[0, 0]:7.4f}) + h{i_hb + 1}@{self.centers[b].name}{b + 1} * ({cf[1, 0]:7.4f})"
                        info2 = f"h{i_ha + 1}@{self.centers[a].name}{a + 1} * ({cf[0, 1]:7.4f}) + h{i_hb + 1}@{self.centers[b].name}{b + 1} * ({cf[1, 1]:7.4f})"
                        
                        io = np.cos(np.arctan2(2 * d12, np.abs(d11 - d22)))
                        print(f"{i_lpo - 1:5d}\t{f'(BD) {host_atoms}, Io = {io:.4f}':<35s}\t{vals[0]:9.5f}\t{info1}\t")
                        print(f"{i_lpo:5d}\t{' ' + host_atoms + ', antibonding (NB)':<35s}\t{vals[1]:9.5f}\t{info2}\t")
                        l2_sum += vals[0] * vals[0]
                        num_bd += 1
                        
        if ry_print_threshold > 0:
            print()
            if num_ry_below_thresh > 0:
                print(f"Note: {num_ry_below_thresh} one-center RY orbitals, each having occupancy below {ry_print_threshold:.2e},")
                print(f" were not printed (use the program option RyOccPrintThreshold to change this behavior)")
                print(f" Total occupancy of these RY orbitals is {ry_below_thresh_occs:.5f}")
                
        print()
        print()
        if lo_info is not None:
            lo_info.bd_per_atomic_pair = self.num_bonds()
            self.print_bond_matrix(lo_info.bd_per_atomic_pair)
            
        if lo_name_str == "CLPO":
            for a in range(self.n_atoms):
                print(f"VAL: \t {self.centers[a].name:2s} \t {lo_info.bd_per_atomic_pair[a, a]:d} ")
                
        print()
        print(f">> {lo_name_str} occupancy summary >>")
        print(f" bonding (BD): {bd_sum:12.5f} in {num_bd:4d} orbitals")
        print(f" anti-bonding (NB): {nb_sum:12.5f} in {num_bd:4d} orbitals")
        print(f" 1c-lone pairs (LP): {lp_sum:12.5f} in {num_lp:4d} orbitals")
        print(f" 1c-unoccupied (RY): {ry_sum:12.5f} in {num_ry:4d} orbitals")
        print()
        print(f"Method BD+LP....in NB+RY BD+NB+LP BD+NB+LP+RY trace(D) Sum[Bd^2+Lp^2] ||D||^2")
        tmp = np.linalg.norm(d_aho, 'fro')
        print(f" {lo_name_str:5s} {bd_sum + lp_sum:12.5f} {num_bd + num_lp:5d} {nb_sum + ry_sum:10.5f} {bd_sum + lp_sum + nb_sum:12.5f} {bd_sum + lp_sum + nb_sum + ry_sum:12.5f} {np.trace(d_aho):12.3f} {l2_sum:12.4f} {tmp * tmp:12.4f}")
        print()
        return d_aho

    def targ_func1(self, hybr2i_h):
        dnew = self.sds_in_hybrid_basis(hybr2i_h, None)
        tmp = 0.0
        for a in range(self.n_atoms):
            for ha in range(self.hybrids_of_atoms[a].n_valid_hybrids):
                ia = hybr2i_h[a][ha]
                b = self.hybrids_of_atoms[a].friend_atom_index[ha]
                tmp += dnew[ia, ia] * dnew[ia, ia]
                if b != -1:
                    hb = self.hybrids_of_atoms[a].friend_hybrid_index[ha]
                    ib = hybr2i_h[b][hb]
                    tmp += dnew[ia, ib] * dnew[ib, ia]
        return tmp

    def iterative_hybrid_opt(self, hybr2i_h, allow_all, iter_comment):
        win_prev = 0.0
        iter_count = 0
        while iter_count < self.n_naos:
            dnew = self.sds_in_hybrid_basis(hybr2i_h, None)
            win_new = self.reconnect_hybrids(dnew, hybr2i_h, allow_all)
            print(f"{iter_comment} iteration {iter_count + 1:2d}: hybrids reconnected, target function = {win_new:12.7f}")
            
            if (iter_count > 0) and (win_new < win_prev - 1.0e-10):
                print(f"WOW!!! new = {win_new:.12e} < prev = {win_prev:.12e} ; new-prev = {win_new - win_prev:.12e}")
                print(f"PHI = {self.targ_func1(hybr2i_h):.12e} tr(DNew**2) = {np.trace(dnew @ dnew.T):.12e} ")
                
            if (iter_count > 0) and (np.abs(win_prev - win_new) < self.target_function_conv_thresh):
                print("Done! ")
                break
                
            win_prev = self.optimize_hybrids()
            print(f"{iter_comment} iteration {iter_count + 1:2d}: hybrids optimized, target function = {win_prev:12.7f}")
            iter_count += 1
            
        print(f"(in {iter_count} iterations)")
        print()

    def create_hybrid_labels(self, prefix, addr2i_h):
        result = [""] * self.n_naos
        for a in range(self.n_atoms):
            for ha in range(len(self.hybrids_of_atoms[a].nao_indices)):
                i_h = addr2i_h[a][ha]
                result[i_h] = f"{self.centers[a].name}{a + 1}:{prefix}{i_h + 1}"
        return result

    def create_clpos_main(self, sds_nao, naos, centers, options):
        self.options = options
        self.target_function_conv_thresh = options.hybr_opt_conv_thresh
        self.opt_max_iter = options.hybr_opt_max_iter
        self.n_atoms = len(centers)
        self.centers = centers
        self.sds_nao = sds_nao
        self.n_naos = len(naos)
        
        print("Creating LPOs (Localized Property-optimized Orbitals)\n\n")
        self.naos_at_center = self.get_naos_at_centers(naos)
        self.hybrids_of_atoms = [AtomicHybrids(self.naos_at_center[a]) for a in range(self.n_atoms)]
        self.dab = self.diatomic_submatrices()
        
        print("Creating initial guess...")
        self.cs_guess()
        
        hybr2i_h = self.get_global_hybrid_indices()
        print("Optimizing LPOs...\n")
        self.opt_lewis_mode = False
        self.iterative_hybrid_opt(hybr2i_h, True, "LPO")
        
        self.lpo_descript = LODescription(self.n_naos)
        self.create_los("LPO ", self.lpo_descript)
        
        hybr2i_h = self.get_global_hybrid_indices()
        print()
        print("*" * 50)
        print("Creating CLPOs (Chemist's Localized Property-optimized Orbitals)\n\n")
        print("Optimizing CLPOs...")
        self.opt_lewis_mode = True
        self.iterative_hybrid_opt(hybr2i_h, False, "CLPO")
        
        hybr2i_h = self.get_global_hybrid_indices()
        self.clpo_descript = LODescription(self.n_naos)
        self.create_los("CLPO", self.clpo_descript)
        self.clpo_descript.hybrid_labels = self.create_hybrid_labels("LHO", hybr2i_h)

    def lo_connectivity(self, los, atomic_charges):
        fragment_of_atom = np.zeros(self.n_atoms, dtype=int)
        valences = np.zeros(self.n_atoms, dtype=int)
        neighs = [np.zeros(self.n_atoms, dtype=int) for _ in range(self.n_atoms)]
        n_neighs = np.zeros(self.n_atoms, dtype=int)
        
        for a in range(self.n_atoms):
            for b in range(a + 1, self.n_atoms):
                if los.bd_per_atomic_pair[a, b] > 0:
                    valences[a] += 1
                    valences[b] += 1
                    neighs[a][n_neighs[a]] = b
                    n_neighs[a] += 1
                    neighs[b][n_neighs[b]] = a
                    n_neighs[b] += 1
                    
        yet_to_visit = np.zeros(self.n_atoms, dtype=int)
        len_yet_to_visit = 0
        fragment_of_atom[:] = -1
        current_fragment_id = 0
        
        for i in range(self.n_atoms):
            if fragment_of_atom[i] != -1: continue
            yet_to_visit[len_yet_to_visit] = i
            len_yet_to_visit += 1
            current_fragment_id += 1
            fragment_of_atom[i] = current_fragment_id
            
            while len_yet_to_visit > 0:
                len_yet_to_visit -= 1
                c = yet_to_visit[len_yet_to_visit]
                for j in range(n_neighs[c]):
                    a2 = neighs[c][j]
                    if fragment_of_atom[a2] == -1:
                        yet_to_visit[len_yet_to_visit] = a2
                        len_yet_to_visit += 1
                        fragment_of_atom[a2] = current_fragment_id
                        
        print(f"There are {current_fragment_id} molecule(s) in the system")
        print("(the 'molecule' is defined as the set of atoms linked with BD orbitals)")
        
        fragment_charges = np.zeros(current_fragment_id)
        n_atoms_in_mol = np.zeros(current_fragment_id, dtype=int)
        for i in range(self.n_atoms):
            f = fragment_of_atom[i] - 1
            fragment_charges[f] += atomic_charges[i]
            n_atoms_in_mol[f] += 1
            
        print()
        print("Molecule TotalNPA NumberOf ListOf ")
        print("id charge atoms atoms...")
        for f in range(current_fragment_id):
            print(f"{f + 1:4d} {fragment_charges[f]:+11.5f} {n_atoms_in_mol[f]:7d} ", end="")
            for i in range(self.n_atoms):
                if fragment_of_atom[i] == f + 1:
                    print(f" {self.centers[i].name}{i + 1}", end="")
            print()
        print()
        print("Note: The total NPA charge is computed as the sum of NPA charge of atoms belonging to each molecule")
        return fragment_of_atom


def create_clpos(
    sds_nao: np.ndarray,
    naos: list[BasisFunction],
    centers: list[AtomicCenter],
    options: CLPOOptions | None = None,
) -> CLPOResult:
    """Full LPO → CLPO pipeline. Replaces PropertyOptimizedOrbitals.createCLPOs(). """
    if options is None:
        options = CLPOOptions()
    engine = PropertyOptimizedOrbitals()
    engine.create_clpos_main(sds_nao, naos, centers, options)
    return CLPOResult(
        lpo=engine.lpo_descript,
        clpo=engine.clpo_descript,
        sds_in_hybrid_basis=engine.sds_in_hybrid_basis(engine.get_global_hybrid_indices(), None)
    )