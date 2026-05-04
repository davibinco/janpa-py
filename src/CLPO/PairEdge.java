package CLPO;

public class PairEdge {
    int a;
    int b;
    double wt;
    PairEdge(int a, int b, double wt) {
        this.a = a;
        this.b = b;
        this.wt = wt;
    }
    @Override
    public String toString() {
        return String.format("PairEdge(a=%d, b=%d, wt=%.2f)", a, b, wt);
    }
}
