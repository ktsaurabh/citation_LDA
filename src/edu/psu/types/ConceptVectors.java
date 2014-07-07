package edu.psu.types;

import java.util.Vector;
import gnu.trove.TIntDoubleHashMap;

public class ConceptVectors extends Document {
	public TIntDoubleHashMap content;
	public Vector<TIntDoubleHashMap> contexts;
	public int docid;
	

	public void setContent(int docid, String name,TIntDoubleHashMap content) {
		this.content = content;
		this.docid = docid;
		this.docName = name;
	}

	public void addContext(TIntDoubleHashMap context) {
		if (this.contexts == null) {
			this.contexts = new Vector<TIntDoubleHashMap>();
		}
		this.contexts.add(context);
	}

	public void processDocuemnt(ContextDocument doc) {
		// process all the sentences
		int radius = 3;
		int[] words = doc.wordSequence.getFeatures();
		for (int i = 0; i < doc.sentenseLength; i++) {
			if (doc.contextObject.contains(i)) {
				TIntDoubleHashMap cxt = new TIntDoubleHashMap();
				int begin = i - radius, end = i + radius;
				// generate context vector
				if (i < radius) {
					begin = i;
				}
				if (i + radius >= doc.sentenseLength) {
					end = doc.sentenseLength - 1;
				}
				for (int j = begin; j <= end; j++) {
					for (int k = doc.SentenseBoundaries.get(j).firstElement(); k < doc.SentenseBoundaries
							.get(j).lastElement(); k++) {
						cxt.adjustOrPutValue(words[k], 1, 1);

					}
				}

				double norm = 0;

				int[] keys = cxt.keys();
				for (int j = 0; j < keys.length; j++) {
					norm += Math.pow(cxt.get(keys[j]), 2);
				}
				norm = Math.pow(norm, 0.5);
				for (int j = 0; j < keys.length; j++) {
					cxt.put(keys[j], cxt.get(keys[j]) / norm);
				}

				for (int citations = 0; citations < doc.contextObject.get(i)
						.size(); citations++) {
					int citationid = doc.contextObject.get(i).get(citations);

					if (!Corpus.concepts.containsKey(citationid)) {
						ConceptVectors cv = new ConceptVectors();
						cv.addContext(cxt);
						Corpus.concepts.put(citationid, cv);
					} else {
						Corpus.concepts.get(citationid).addContext(cxt);
					}
				}
			}
		}

	}

}