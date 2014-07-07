package edu.psu.types;

import java.util.*;

import gnu.trove.*;
import edu.psu.topic.Model;
import edu.psu.types.*;
import java.io.*;

public class Document implements Serializable {

	public LabelSequence topicAssignments = null; // sequence representation of
													// topic assignments
	public FeatureSequence wordSequence = null; // sequence representation of
												// words in the document
	public TIntIntHashMap topic_counts = null; // bag representation of topic
												// assignments
	public TIntIntHashMap word_counts = null; // bag representation of words in
												// document
	public TIntDoubleHashMap tfidfVector = null; // bag representation of words
													// in
	// document
	public TIntDoubleHashMap normed_tfidfVector = null;

	public TIntDoubleHashMap normed_tfVector = null;
	public TIntIntHashMap links = null;
	public Vector<Integer> linksTopics = null;
	public String docName;
	public int docId;
	public int docLabel;
	public int docLength;
	public double norm_tfidf;
	public double norm_tf;

	public Document() {
	}
	
	public Document(String name, String line) {
		String[] flds = line.split(":");
		if(flds.length==2){
			String[] content = flds[0].split(" ");
			String[] citations = flds[1].split(" ");
			ProcessCSV(name, 0, content,citations);
		}
		if(flds.length==1){
			String[] content = flds[0].split(" ");
			String[] citations = null;
			ProcessCSV(name, 0, content,citations);
		}
	}

	public Document(LabelSequence ts, FeatureSequence fs) {
		// assuming the content has been processed
		topicAssignments = ts;
		wordSequence = fs;
	}
	
	
	public void addDocument(String content, String links) {
		// create document object from content
		String name;
		String[] words, citations = null;
		int label;
		String[] contents = content.split(" ");
		name = contents[0];
		label = Integer.parseInt(contents[1]);
		words = new String[contents.length - 2];
		for (int i = 2; i < contents.length; i++) {
			words[i - 2] = contents[i];
		}
		
		//To-do add citations
		ProcessCSV(name, label, words, citations);

	}

	public void ProcessCSV(String name, int label, String[] content,
			String[] citations) {

		// create document object from content
		docName = name;
		docId = Corpus.size;
		docLabel = label;
		docLength = content.length;
		wordSequence = new FeatureSequence(Corpus.vocabulary);
		word_counts = new TIntIntHashMap();
		links = new TIntIntHashMap();

		// String[] words = content.split(" ");
		// String[] edges = citations.split(" ");
		if (content == null) {
			// System.out.println("INFO: No content");

		} else {
			for (int i = 0; i < content.length; i++) {
				String word = content[i].trim();
				int index = Corpus.vocabulary.lookupIndex(word, true);

				wordSequence.add(index);

			}
			int[] words = wordSequence.getFeatures();
			for (int i = 0; i < wordSequence.length; i++) {
				word_counts.adjustOrPutValue(words[i], 1, 1);
			}
		}
		if (citations == null) {
			// System.out.println("INFO: No citations");

		} else {
			for (int i = 0; i < citations.length; i++) {
				String edge = citations[i].trim();
				int index = Corpus.docAlphabet.lookupIndex(edge, true);
				//links.add(index);
				links.adjustOrPutValue(Corpus.docAlphabet.lookupIndex(edge), 1, 1);

			}
		}
	}
	
	public void ProcessCSV(String name, int label, int docid, Vector<String> content,
			Vector<String> citations) {

		// create document object from content
		docName = name;
		docId = docid;
		docLabel = label;
		docLength = content.size();
		wordSequence = new FeatureSequence(Corpus.vocabulary);
		word_counts = new TIntIntHashMap();
		links = new TIntIntHashMap();

		// String[] words = content.split(" ");
		// String[] edges = citations.split(" ");
		if (content == null) {
			// System.out.println("INFO: No content");

		} else {
			for (int i = 0; i < content.size(); i++) {
				String word = content.get(i);
				int index = Corpus.vocabulary.lookupIndex(word, true);

				wordSequence.add(index);

			}
			int[] words = wordSequence.getFeatures();
			for (int i = 0; i < wordSequence.length; i++) {
				word_counts.adjustOrPutValue(words[i], 1, 1);
			}
		}
		if (citations == null) {
			// System.out.println("INFO: No citations");

		} else {
			for (int i = 0; i < citations.size(); i++) {
				String edge = citations.get(i);
				int index = Corpus.docAlphabet.lookupIndex(edge);
				//links.add(index);
				links.adjustOrPutValue(Corpus.docAlphabet.lookupIndex(edge), 1, 1);

			}
		}
	}

	public void normalizeTFIDF() {
		if (wordSequence == null || wordSequence.length == 0) {
			System.out.println("WordSequence Empty for Document " + docName);
			return;
		}
		tfidfVector = new TIntDoubleHashMap();
		normed_tfidfVector = new TIntDoubleHashMap();
		normed_tfVector = new TIntDoubleHashMap();

		int[] words = word_counts.keys();
		for (int i = 0; i < words.length; i++) {
			tfidfVector.put(words[i],
					word_counts.get(words[i]) * Corpus.idf.get(words[i]));

		}
		norm_tfidf = 0.0;
		for (int i = 0; i < words.length; i++) {
			norm_tfidf += Math.pow(tfidfVector.get(words[i]), 2.0);
		}
		norm_tfidf = Math.pow(norm_tfidf, 0.5);
		for (int i = 0; i < words.length; i++) {
			normed_tfidfVector.put(words[i], tfidfVector.get(words[i])
					/ norm_tfidf);
		}
		norm_tf = 0.0;
		for (int i = 0; i < Corpus.vocabulary.size(); i++) {
		//for (int i = 0; i < words.length; i++) {
			norm_tf += Math.pow((word_counts.get(i)-Corpus.mean[i]), 2.0);
			//norm_tf += Math.pow((word_counts.get(i)), 2.0);

		}
		norm_tf = Math.pow(norm_tf, 0.5);
		 for (int i = 0; i < Corpus.vocabulary.size(); i++) {
		//for (int i = 0; i < words.length; i++) {
			normed_tfVector.put(i, (word_counts.get(i) - Corpus.mean[i])/
			 norm_tf);

		//	normed_tfVector.put(i, ((double)word_counts.get(i)) / norm_tf);
		}

	}

	private static final long serialVersionUID = 1;
	private static final int CURRENT_SERIAL_VERSION = 0;

	private void writeObject(ObjectOutputStream out) throws IOException {

		out.writeInt(CURRENT_SERIAL_VERSION);

		out.writeInt(docLabel);
		out.writeInt(docId);
		out.writeInt(docLength);
		out.writeDouble(norm_tfidf);
		out.writeObject(docName);
		out.writeObject(topicAssignments);
		out.writeObject(wordSequence);
		out.writeObject(tfidfVector);
		out.writeObject(normed_tfidfVector);
		out.writeObject(topic_counts);
		out.writeObject(word_counts);
		out.writeObject(links);
		out.writeObject(linksTopics);

	}

	private void readObject(ObjectInputStream in) throws IOException,
			ClassNotFoundException {
		int version = in.readInt();

		docLabel = in.readInt();

		docId = in.readInt();

		docLength = in.readInt();

		norm_tfidf = in.readDouble();

		docName = in.readObject().toString();

		topicAssignments = (LabelSequence) in.readObject();
		wordSequence = (FeatureSequence) in.readObject();
		tfidfVector = (TIntDoubleHashMap) in.readObject();
		normed_tfidfVector = (TIntDoubleHashMap) in.readObject();
		topic_counts = (TIntIntHashMap) in.readObject();
		word_counts = (TIntIntHashMap) in.readObject();
		links = (TIntIntHashMap) in.readObject();
		linksTopics = (Vector<Integer>) in.readObject();

	}

}


