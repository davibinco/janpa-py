package ProgramOptions;

import java.io.*;
import ProgramOptions.WarningManager;

/**
 * A class for convenient automatic parameter / command line parsing,
 * storing program work options and accessing them by name / by object
 * 
 * This file is a part of the JANPA project. 
 *
 * Copyright (c) 2014, Tymofii Nikolaienko
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 
 * 3. All advertising and/or published materials mentioning features or use 
 *    of this software must display the following acknowledgement:
 * 
 *       This product includes components from JANPA package of programs
 *       ( http://janpa.sourceforge.net/ ) developed by Tymofii Nikolaienko
 * 
 * 4. Neither the name of the developer, Tymofii Nikolaienko,  nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 * 
 * 5. In case if the code of JANPA package code and/or its parts and/or any data 
 *    produced with JANPA package of programs are published, the following citations
 *    for the JANPA package of programs should be given:
 *
 *     1) T.Y.Nikolaienko, L.A.Bulavin; Int. J. Quantum Chem. (2019), 
 *        Vol.119, page e25798, DOI: 10.1002/qua.25798
 *     2) T.Y.Nikolaienko, L.A.Bulavin, D.M.Hovorun; Comput.Theor.Chem. (2014),
 *        V. 1050, P. 15-22, DOI: 10.1016/j.comptc.2014.10.002
 * 
 * THIS SOFTWARE IS PROVIDED BY ''AS IS'' AND ANY EXPRESS OR IMPLIED 
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
 * IN NO EVENT SHALL Tymofii Nikolaienko BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 *
 * @author (c) Tymofii Nikolaienko, 2014
 */

/* TODO: add support for object classes String[], Double[], ...
  (for arguments like "-arr 1.0 -arr 2.0 -arr 3.0" => double[]{1.0,2.0,3.0})
 * 
 */

public class OptionParameter {
    // Some (private) fields
    private OptionParameter _next = null;
    private Object data;
    public String Name;
    public String Description;
    //----------------------------------------------------------------------
    // default formats for data->text conversion
    public static String default_double_format = " %.20E";
    public static String default_int_format = "%d";
    //----------------------------------------------------------------------
    // A stream for printing messages / warnings
    public PrintStream out = System.out;
    //----------------------------------------------------------------------
    /** Initializes local fields and is called from constructors */
    private void _init(OptionParameter FirstInList, String Name,Object data, String Description) {
        if (FirstInList == null)
            this._next = null; // we're beginning the list
        else {
            // append to list
            OptionParameter tmp = FirstInList;
            while (tmp._next != null) tmp = tmp._next;
            tmp._next = this;
        }
        this.Name = Name;
        this.data = data;
        this.Description = Description;
    }
    //----------------------------------------------------------------------
    /** A constructor:
     * Creates an object AND appends it to the list beginning with FirstInList
     * in case if FirstInList != null
     * 
     * @param data: a data object String / Double / Boolean / ...
     *  It is good idea NOT to pass null as the data
     */
    public OptionParameter(OptionParameter FirstInList, String Name,Object data, String Description) {
        _init(FirstInList, Name, data, Description);
    }
    // and its simplieifed version
    public OptionParameter(OptionParameter FirstInList, String Name,Object data) {
        _init(FirstInList, Name, data, "");
    }
    //----------------------------------------------------------------------
    void On_Error(String message) {
        System.out.println("ERROR " + message);
    }
    //----------------------------------------------------------------------
    // Methods for retrieving value
    //----------------------------------------------------------------------
    /** returns 'raw' object (not recommended for use other than instanceof-checking) */
    public Object get_data() {
        return data;
    }
    //----------------------------------------------------------------------
    public String get_String() {//throws Exception {
        if (!(data instanceof String)) On_Error(String.format("get_String() is impossible for option \"%s\" with Class.getName() = %s", Name, data.getClass().getName()));
            //throw new Exception();
        return (String)data;
    }
    //----------------------------------------------------------------------
    public double get_double() {//throws Exception {
        if (!(data instanceof Double)) On_Error(String.format("get_double() is impossible for option \"%s\" with Class.getName() = %s", Name, data.getClass().getName()));
            //throw new Exception(String.format("get_double() is impossible for option \"%s\" with Class.getName() = %s", Name, data.getClass().getName()));
        return ((Double)data).doubleValue();
    }
    //----------------------------------------------------------------------
    public boolean get_boolean() {//throws Exception {
        if (!(data instanceof Boolean)) On_Error(String.format("get_boolean() is impossible for option \"%s\" with Class.getName() = %s", Name, data.getClass().getName()));
            //throw new Exception(String.format("get_boolean() is impossible for option \"%s\" with Class.getName() = %s", Name, data.getClass().getName()));
        return ((Boolean)data).booleanValue();
    }
    /*
    public boolean get_boolean_noExceptions() {
        return ((Boolean)data).booleanValue(); // if data is not isinstanceof, runtime exception will stop the program
    }
     * 
     */
    //----------------------------------------------------------------------
    public int get_int() {//throws Exception {
        if (!(data instanceof Integer)) On_Error(String.format("get_int() is impossible for option \"%s\" with Class.getName() = %s", Name, data.getClass().getName()));
            //throw new Exception(String.format("get_int() is impossible for option \"%s\" with Class.getName() = %s", Name, data.getClass().getName()));
        return ((Integer)data).intValue();
    }
    //----------------------------------------------------------------------
    /** @returns A human-readable representation of any known kind of data
     * An internal version
     */
    public static String get_data_as_String(Object data) throws Exception {
            if (data instanceof String) {
                if (((String)data).isEmpty())       // upd. 10.01.2016
                    return "(empty)";
                else
                    return (String)data; 
            } else
            if (data instanceof Boolean) return ((Boolean)data).toString(); else
            if (data instanceof Double) return String.format(default_double_format, (Double)data); else // using get_double() would be better, but get_double() is not static
            if (data instanceof Integer) return String.format(default_int_format, (Integer)data); else // the same with get_int()
            if (data instanceof Object[]) {
                // an array of objects of any known type
                String result = "{ ";
                Object[] A = (Object[])data;
                for (int i=0; i<A.length; ++i) {
                    if (i>0) result += ", ";
                    result += get_data_as_String(A[i]);
                }
                return result + " }";
            }
                //On_Error(" The class of data is unknown!");
                throw new Exception(" The class of data is unknown!");
    }
    //----------------------------------------------------------------------
    /**
     * A public version of get_data_as_String(), which converts any type to String
     */
    public String get_as_String() /*throws Exception*/ { 
        try {
            return get_data_as_String(data);
        } catch (Exception e) {
            On_Error(e.getMessage());
            return null;
        }
    }
    //----------------------------------------------------------------------
    /**
     * Prints pairs (parameter)-(value) for the list beginning with @param FirstInList
     * using a specified printf-compatible @param format
     * If FirstInList == null, this is assumed to be FirstInList
     * TODO: for arrays the output of this method can not be parsed by ReadDataFromString
     */
    public void Print_as_List(OptionParameter FirstInList, String format) throws Exception {
        if (FirstInList == null) FirstInList = this;        // Note: Java passes everything by value only!
        OptionParameter tmp = FirstInList;
        if (format == null) format = "%s   %s%n";
        while (tmp != null) {
            if (!((Object)tmp instanceof Object[]))
                // for non-array types
                out.printf(format, tmp.Name, tmp.get_as_String());
            else {
                // for array types
                Object[] A = (Object[])data;
                for (int i=0; i<A.length; ++i) {
                    out.printf(format, tmp.Name, get_data_as_String(A[i])); // print the same name near each element of the array
                    // obviously, this simple loop will not print nested arrays appropriately // TODO
                }
            }
            tmp = tmp._next;
        }
        out.flush();
    }
    //----------------------------------------------------------------------
    /**
     * Prints a list of (parameter name)-(description) (one pair per line)
     * for each member of the list beginning with @param FirstInList;
     * formats each pair using a @param format parameter (if it is null,
     * internal template is used).
     * 
     */
    public void PrintDescriptions(OptionParameter FirstInList, String format) {//throws Exception{
        if (FirstInList == null) FirstInList = this;        // Note: Java passes everything by value only!
        if (format == null) format = "%s  %s%n"; // up to 3 %s are available
        while (FirstInList != null) {
            out.printf(format, FirstInList.Name, FirstInList.Description, FirstInList.get_as_String());
            FirstInList = FirstInList._next;
        }
        
    }
    //----------------------------------------------------------------------

    // Methods for changing data
    //----------------------------------------------------------------------
    /** Sets 'raw' data.
     * This method should be used like (...).set_data(new Double())
     */
    public void set_data(Object data) {
        this.data = data;
    }
    //----------------------------------------------------------------------
    /** An internal parsing routine used by the public method ReadDataFromString()
     @returns the Object of the same class as old_data with the proper value read from the String s
     * sets datatypoe_known[0] to true if the data type (passed in old_data) is known; 
     * Note that datatypoe_known is declared as array only in order to make datatypoe_known[0] to 
     * be changeable from the method (Java passes everything by value!)
     *
     */
    private Object _parse_data(String s, Object old_data, boolean[] datatype_known /*holds a data type and gets a result*/) {
        datatype_known[0] = true;
        if (old_data instanceof Boolean)
            // the presence of the key with a given name is sufficient for setting the value to TRUE
            return (s == null)?(true):(Boolean.valueOf(s)); //the same as 'if (s == null) data = true; else data = Boolean.valueOf(s);'
        else
        if (old_data instanceof String)  return (s==null)?(""):(s /*no 'new String()' is needed for Java*/); else
        if (old_data instanceof Double)  return (s==null)?(new Double(0.0)):(new Double(s)); else
        if (old_data instanceof Integer) return (s==null)?(new Integer(0)):(new Integer(s)); else
            datatype_known[0] = false;
        return old_data;
    }
    //----------------------------------------------------------------------
    /**
     *      A KEY PUBLIC METHOD FOR READING OBJECTS FROM TEXT DATA
     * Parses string @param s according to the parameter type
     * if s = null, defauls values are used (String: "", Integer/Double: 0, Boolean: true)
     */
    public void ReadDataFromString(String s) throws Exception {
        boolean[] datatype_known = new boolean[]{true};
        // invoke a parser for non-array types
        data = _parse_data(s, data, datatype_known);
        if (!datatype_known[0]) {
            // if not a simple, 1-element data type => array or unknown
            if (data instanceof Object[]) {
                // a special case for treating arrays
                if (s == null)
                    //data = new Object[0]; // clear the array - it is not possible since the datatype will be lost!
                    out.println(" ERROR: value can not be empty for "+this.Name);
                else {
                    // add an element to the array
                    Object[] A = (Object[])data;
                    Object[] new_data = new Object[A.length+1];
                    System.arraycopy((Object[])data, 0, new_data, 0, A.length);
                    new_data[A.length] = _parse_data( s, A[A.length-1] /*all elements of the array should have the same type!*/, datatype_known );
                }
            } else
                throw new Exception(" The class of data is unknown!");
        }
        
    }

    //----------------------------------------------------------------------
    /**
     * Returns an object for a parameter with a given name if it is found in the list
     * beginning with FirstInList;
     * Returns null if none found.
     * If FirstInList == null, this is assumed to be FirstInList
     */
    public OptionParameter FindInListByName(OptionParameter FirstInList, String name, boolean ignore_case)  {
        if (FirstInList == null) FirstInList = this;        // Note: Java passes everything by value only!
        OptionParameter result = FirstInList;
        boolean found = false;
        while ((result != null) && (!found)) {
            if (ignore_case)
                found = name.toLowerCase().equals(result.Name.toLowerCase());
            else
                found = name.equals(result.Name);

            if (!found)
                result = result._next;
        }
        // if (!found) then the result is guaranteed to be null !
        return result;
    }
    //----------------------------------------------------------------------
    /**
     * Parses a command line splitted into an @param args array:
     * Example: args={"name1","value1","name2","name3","value3"}
     * Features:
     * + Once a name of a parameter of type Boolean is encountered,
     *   its value is set to true (i.e., no "value" parameter is required!)
     *    => Booleans can only be set to 'true'
     * + String parameters can not be set to "", since once the name of a String-
     *   type parameter is encountered, the next item of array is assumed
     *   to be a value.
     * + Each name is being searched in the entire list of OptionParameters
     *   starting from FirstInList; if FirstInList == null, it is set to this
     * 
     */
    public void LoadOptionsFromCommandLine(String[] args, OptionParameter FirstInList, boolean ignore_case,
            boolean fail_in_unknown) throws Exception
    {
        //boolean result = true;
        if (FirstInList == null) FirstInList = this;
        int i = 0;
        OptionParameter tmp;
        while (i < args.length) {
            tmp = FirstInList.FindInListByName(null, args[i], ignore_case);
            if (tmp == null) {
                if (fail_in_unknown) 
                    throw new Exception("parameter \"" + args[i] + "\" is unknown!%n");
                WarningManager.warning_printf("Warning: parameter \"%s\" is unknown!%n",args[i]);
                ++i; // go to the next enrty
            } else {
                // if the parameter is known
                if (tmp.get_data() instanceof Boolean) {
                    tmp.set_data(true);
                    ++i;
                } else {
                    // all other types except of Boolean require the next entry to
                    // be the parameter value
                    ++i; // go to the next enrty
                    if (i<args.length) {
                        tmp.ReadDataFromString(args[i]);
                        ++i; // go to the next enrty
                    } else
                        out.printf("ERROR: value missing for parameter \"%s\"; ignoring%n",args[i-1]);
                }
            }
        }
        //return result;
    }
    //----------------------------------------------------------------------
    /**
     * Works similar to the LoadOptionsFromCommandLine, but each line is limited
     * to only one Name-Value pair.
     * Therefore, 
     *  + booleans can be set to false (if a corresponding line contains more
     *    than just a parameter name, but the "false" value also);
     *  + Strings can be set to ""
     * Features:
     *  + One or more spaces (" ") delimit parameter name from its value; all spaces at the
     *    line beginning are ignored; 
     *  + All lines which do not begin with a known parameter name are 'quietly' ignored
     *  + Entire file is being read
     */
    public void LoadOptionsFromFile(BufferedReader fr, OptionParameter FirstInList, boolean ignore_case) throws Exception {
        if (FirstInList == null) FirstInList = this;
        String s;
        final char delimiter = ' '; // a signal of the 'name' field end
        while ((s = fr.readLine()) != null) {
            // skip the spaces/tabs at the line beginning
            s = s.replaceFirst("[ \\t]*", "");
            if (!s.isEmpty()) {
                int pos1 = s.indexOf(delimiter); // first space after the parameter name
                String name = null, value = null;
                if (pos1 == (-1)) 
                    // a space has not been found => s = "name"
                    name = s;
                else {
                    // a space HAS been found => s = "name " or s = "name value"
                    name = s.substring(0, pos1);
                    if (pos1 < s.length()) {
                        // s = "name value"
                        value = s.substring(pos1+1).replaceFirst("[ \\t]*", ""); // skip all spaces/tabs at the beginning of the value
                        if (value.isEmpty()) value = null; // since, e.g., boolean with the value "" will not be treated correctly
                    }
                    // otherwise, value = null
                }
                OptionParameter tmp = FirstInList.FindInListByName(null, name, ignore_case);
                if (tmp != null)
                    tmp.ReadDataFromString(value); // value == null IS acceptable
            }
        }
    }
    //----------------------------------------------------------------------
    /*
    // 1-st index: number, 2-nd index: 0=Name, 1=Value
    // returns true if teh specified name has been found
    public boolean FindAndParse(String[][] name_value_pairs, boolean ignore_case) throws Exception {
        boolean found = false;
        int i=0;
        String lc_Name = Name.toLowerCase();
        // is there our name in the list ?
        while ((!found) && (i<name_value_pairs.length)) {
            if (ignore_case)
                found = name_value_pairs[i][0].toLowerCase().equals(lc_Name);
            else
                found = name_value_pairs[i][0].equals(Name);

            if (!found) ++i;
        }
        // parse the value if the specified name has been found
        if (found) {
            String s = name_value_pairs[i][1];
            this.ReadDataFromString(s);
        }

        return found;
    }
     * 
     */
    //----------------------------------------------------------------------
    }
