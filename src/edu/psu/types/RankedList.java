package edu.psu.types;

import java.util.Iterator;
import java.util.TreeSet;
import java.util.Vector;

import gnu.trove.*;

public class RankedList {
	public Vector<Integer> indices;
	public TIntDoubleHashMap values;

	public RankedList() {
		indices = new Vector<Integer>();
		values = new TIntDoubleHashMap();
	}

	public void add(int id, double value) {
		values.adjustOrPutValue(id, value, value);
	}

	public void sortListDesc() {
		int[] keys = values.keys();
		TreeSet<IDSorter> sortedCitation = new TreeSet<IDSorter>();
		for (int i = 0; i < keys.length; i++) {
			sortedCitation.add(new IDSorter(keys[i], values.get(keys[i])));
		}
		Iterator<IDSorter> iter = sortedCitation.iterator();
		int index = 0;
		while (iter.hasNext()) {
			IDSorter ids = iter.next();
			indices.add(ids.getID());
		}
		if (indices.size() != values.size()) {
			System.out
					.println("Error in Sorting the list!! size does not match!!");
			System.exit(0);
		}
	}

	public double KindallTauCoefficient(RankedList list) {
		double coeff = 0;
		if (this.indices.size() != list.indices.size()) {
			System.out.println("Size of List 1 = " + this.indices.size()
					+ "; Size of List 2 = " + list.indices.size());
			System.out.println("Lists are not of equal sizes!! ");
			// System.exit(0);
		}
		int conc = 0, disc = 0;
		int[] keys = values.keys();
		TDoubleIntHashMap group1 = new TDoubleIntHashMap();
		for (int i = 0; i < keys.length; i++) {
			if (list.values.containsKey(keys[i])) {
				group1.adjustOrPutValue(values.get(keys[i]), 1, 1);
			}
		}
		keys = list.values.keys();
		TDoubleIntHashMap group2 = new TDoubleIntHashMap();
		for (int i = 0; i < keys.length; i++) {
			if (values.containsKey(keys[i])) {
				group2.adjustOrPutValue(list.values.get(keys[i]),
						1, 1);
			}
		}
		int ties1 = 0, ties2 = 0;
		double[] keys_t = group1.keys();
		for (int i = 0; i < keys_t.length; i++) {
			ties1 += (group1.get(keys_t[i]) * (group1.get(keys_t[i]) - 1) / 2);
		}
		keys_t = group2.keys();
		for (int i = 0; i < keys_t.length; i++) {
			ties2 += (group2.get(keys_t[i]) * (group2.get(keys_t[i]) - 1) / 2);
		}
		double total = 0;
		keys = values.keys();
		for (int i = 0; i < keys.length; i++) {
			if (list.values.containsKey(keys[i])) {
				for (int j = i + 1; j < keys.length; j++) {
					if (list.values.containsKey(keys[j])) {
						total++;
						if ((values.get(keys[i]) != values.get(keys[j]))
								&& list.values.get(keys[i]) != list.values
										.get(keys[j])) {
							if (((values.get(keys[i]) > values.get(keys[j])) && (list.values
									.get(keys[i]) > list.values.get(keys[j])))
									|| ((values.get(keys[i]) < values
											.get(keys[j])) && (list.values
											.get(keys[i]) < list.values
											.get(keys[j])))) {
								conc++;
							} else if (((values.get(keys[i]) > values
									.get(keys[j])) && (list.values.get(keys[i]) < list.values
									.get(keys[j])))
									|| ((values.get(keys[i]) < values
											.get(keys[j])) && (list.values
											.get(keys[i]) > list.values
											.get(keys[j])))) {
								disc++;
							}
						}

					}

				}
			}

		}

		coeff = ((conc - disc) / Math.pow((total - ties1) * (total - ties2),
				0.5));
		return coeff;
	}

}
