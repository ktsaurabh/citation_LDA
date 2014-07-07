package edu.psu.types;


public class IDSorter implements Comparable {
    int id; double p;
    public IDSorter (int id, double p) { this.id = id; this.p = p; }
    public IDSorter (int id, int p) { this.id = id; this.p = p; }
    public final int compareTo (Object o2) {
	if (p > ((IDSorter) o2).p)
	    return -1;
	//else if (p == ((IDSorter) o2).p)
	//    return 0;
	else return 1;
    }
    public int getID() {return id;}
    public double getWeight() {return p;}
    /** Reinitialize an IDSorter */
    public void set(int id, double p) { this.id = id; this.p = p; }
}