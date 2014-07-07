package edu.psu.types;


public class DescSorter implements Comparable {
    String id; double p;
    public DescSorter (String id, double p) { this.id = id; this.p = p; }
    public DescSorter (String id, int p) { this.id = id; this.p = p; }
    public final int compareTo (Object o2) {
	if (p > ((DescSorter) o2).p)
	    return -1;
	else //if (p < ((DescSorter) o2).p)
	    return 1;
	//else return 0;
    }
    public String getID() {return id;}
    public double getWeight() {return p;}
    /** Reinitialize an IDSorter */
    public void set(String id, double p) { this.id = id; this.p = p; }
}