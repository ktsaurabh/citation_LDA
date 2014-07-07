package edu.psu.types;

import gnu.trove.TIntDoubleHashMap;
import edu.psu.types.Corpus;
public class WordOccurrences{
	public int id;
	public TIntDoubleHashMap occurrences;
	public WordOccurrences(int id, Corpus corpus){
		this.id = id;
		occurrences = new TIntDoubleHashMap();
		for(Document doc: corpus.docs){
			if(doc.word_counts.containsKey(id)){
				occurrences.put(doc.docId, doc.normed_tfVector.get(id));
			}
		}
	}
}