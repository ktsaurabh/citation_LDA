package edu.psu.types;

import java.util.*;
import java.util.regex.Pattern;

import gnu.trove.*;
import edu.psu.types.*;
import java.io.*;

/*
 * Data structure
 * Initilization
 * Sampling
 * Results
 * 
 * 
 * 
 * 
 */

public class ContextDocument extends Document {

	// Call docuements constructor to provide content and citations
	class CitationContext {
		public int citationid;
		public int begin;
		public int end;
		public int length;

		public CitationContext(int citationid, int begin, int end) {
			this.citationid = citationid;
			this.begin = begin;
			this.end = end;
			this.length = end - begin + 1;
		}

		public void updateLength() {
			this.length = end - begin + 1;
		}
	}

	// in its constructor
	public TIntObjectHashMap<Vector<Integer>> contextObject = null;
	public TIntObjectHashMap<Vector<Integer>> SentenseBoundaries = null;
	public TIntIntHashMap citationSet;
	public TIntIntHashMap citationSetTest;
	public TIntIntHashMap citationTopicAssignment;
	public int sentenseLength;
	public int wordLength;
	public Vector<CitationContext> citationContexts;
	public double avgContextLength;
	public double avgLeftContextLength;
	public double avgRightContextLength;

	public ContextDocument() {
		wordLength = 0;
		sentenseLength = 0;
		SentenseBoundaries = new TIntObjectHashMap<Vector<Integer>>();
		contextObject = new TIntObjectHashMap<Vector<Integer>>();
		citationSet = new TIntIntHashMap();
		citationSetTest = new TIntIntHashMap();
		citationContexts = new Vector<CitationContext>();

	}

	/*
	 * public ContextDocument(String context){ //Format: sentense:citations\n
	 * String[] flds = context.split(":"); addDocument(flds[0], flds[1]);
	 * 
	 * }
	 */
	public TIntDoubleHashMap normalizeSparseVector(TIntDoubleHashMap doc) {
		double norm = 0;

		int[] keys = doc.keys();
		for (int j = 0; j < keys.length; j++) {
			norm += Math.pow(doc.get(keys[j]), 2);
		}
		norm = Math.pow(norm, 0.5);
		for (int j = 0; j < keys.length; j++) {
			doc.put(keys[j], doc.get(keys[j]) / norm);
		}
		return doc;

	}

	public double topicVectorSimilarity(TIntDoubleHashMap v1,
			TIntDoubleHashMap v2) {
		double sim = 0;
		int[] keys = v1.keys();
		for (int i = 0; i < keys.length; i++) {
			sim += (v1.get(keys[i]) * v2.get(keys[i]));
		}

		return sim;
	}

	public void add(String name, String context) {
		int label = 0;

		Vector<Integer> tmpSentence = new Vector<Integer>();
		Vector<Integer> tmpContext = new Vector<Integer>();
		tmpSentence.add(wordLength);
		// Format: sentense:citations\n
		String[] flds = context.split(":");
		// addDocument(flds[0], flds[1]);
		String[] words = flds[0].split(" ");
		// label = Integer.parseInt(words[1]);
		String[] citations = null;
		if (flds.length > 1) {
			citations = new String[flds.length - 1];
			for (int i = 1; i < flds.length; i++) {
				citations[i - 1] = flds[i];
			}
		}
		// add citations first so that citation id are obtained

		// for (int i = 0; i < citations.length; i++)
		// tmpContext.add();
		ProcessCSV(name, label, words, citations);
		if (citations != null) {
			for (int i = 0; i < citations.length; i++) {
				if (contextObject.contains(sentenseLength)) {
					contextObject.get(sentenseLength).add(
							Corpus.docAlphabet.lookupIndex(citations[i]));

				}

				else {
					tmpContext
							.add(Corpus.docAlphabet.lookupIndex(citations[i]));
					contextObject.put(sentenseLength, tmpContext);
				}
				citationSet
						.put(Corpus.docAlphabet.lookupIndex(citations[i]), 1);
			}
		}
		wordLength += words.length;
		sentenseLength++;
		tmpSentence.add(wordLength - tmpSentence.firstElement());
		SentenseBoundaries.put(sentenseLength, tmpSentence);

	}

	public void add(String name, int docid, Vector<String> context) {
		int label = 0;
		Vector<String> contentVector = new Vector<String>();
		Vector<String> citationVector = new Vector<String>();
		for (int line = 0; line < context.size(); line++) {
			if (Pattern.matches("[a-zA-Z0-9]+", context.get(line)) || true) {
				// System.out.println("Here");
				Vector<Integer> tmpSentence = new Vector<Integer>();
				Vector<Integer> tmpContext = new Vector<Integer>();
				tmpSentence.add(wordLength);
				// Format: sentense:citations\n
				String[] flds = context.get(line).split(":");
				// addDocument(flds[0], flds[1]);
				String[] words = flds[0].split(" ");
				// label = Integer.parseInt(words[1]);
				String[] citations = null;
				if (flds.length > 1) {
					citations = new String[flds.length - 1];
					for (int i = 1; i < flds.length; i++) {
						citations[i - 1] = flds[i];
						// System.out.println(flds[i]);
					}
				}
				// add citations first so that citation id are obtained

				// for (int i = 0; i < citations.length; i++)
				// tmpContext.add();

				if (citations != null) {
					for (int i = 0; i < citations.length; i++) {
						CitationContext cc = new CitationContext(
								Corpus.docAlphabet.lookupIndex(citations[i]),
								sentenseLength, sentenseLength);
						citationContexts.add(cc);
						if (contextObject.contains(sentenseLength)) {
							contextObject.get(sentenseLength).add(
									Corpus.docAlphabet
											.lookupIndex(citations[i]));
						} else {
							tmpContext.add(Corpus.docAlphabet
									.lookupIndex(citations[i]));
							contextObject.put(sentenseLength, tmpContext);
						}
						citationSet
								.put(Corpus.docAlphabet
										.lookupIndex(citations[i]), 1);

					}
				}
				wordLength += words.length;
				tmpSentence.add(wordLength);
				SentenseBoundaries.put(sentenseLength, tmpSentence);
				sentenseLength++;
				for (int i = 0; i < words.length; i++) {
					contentVector.add(words[i]);
				}
				if (citations != null) {
					for (int i = 0; i < citations.length; i++) {
						citationVector.add(citations[i]);
					}
				}
			}
		}
		ProcessCSV(name, label, docid, contentVector, citationVector);
		if (wordLength != wordSequence.length) {
			System.out.println("wordSequence.length=" + wordSequence.length);
			System.out.println("wordLength=" + wordLength);
		}
		// add concept vectors
		double norm = 0;
		TIntDoubleHashMap content = new TIntDoubleHashMap();
		int[] keys = this.word_counts.keys();
		for (int i = 0; i < keys.length; i++) {
			norm += Math.pow(this.word_counts.get(keys[i]), 2);
		}
		norm = Math.pow(norm, 0.5);
		for (int i = 0; i < keys.length; i++) {
			content.put(keys[i], this.word_counts.get(keys[i]) / norm);
		}
		//

		if (!Corpus.concepts.containsKey(docid)) {
			ConceptVectors cv = new ConceptVectors();
			cv.setContent(docid, name, content);
			Corpus.concepts.put(docid, cv);
			cv.processDocuemnt(this);
		} else {
			Corpus.concepts.get(docid).setContent(docid, name, content);
			Corpus.concepts.get(docid).processDocuemnt(this);
		}

	}

	public ContextDocument RemoveCitation() {
		ContextDocument doc = new ContextDocument();
		doc.word_counts = this.word_counts;
		doc.wordSequence = this.wordSequence;
		doc.topic_counts = this.topic_counts;
		doc.topicAssignments = this.topicAssignments;
		doc.docId = this.docId;
		doc.docName = this.docName;
		doc.docLength = this.docLength;
		doc.docLabel = this.docLabel;
		doc.wordLength = this.wordLength;
		doc.sentenseLength = this.sentenseLength;
		int[] keys = this.SentenseBoundaries.keys();
		for (int i = 0; i < keys.length; i++) {
			Vector<Integer> sentences = new Vector<Integer>();
			for (int j = 0; j < this.SentenseBoundaries.get(i).size(); j++) {
				sentences.add(this.SentenseBoundaries.get(i).get(j));
			}
			doc.SentenseBoundaries.put(keys[i], sentences);
		}
		// copy the citationset in a seperate variable
		keys = this.citationSet.keys();
		for (int i = 0; i < keys.length; i++) {
			doc.citationSetTest.adjustOrPutValue(keys[i], 1, 1);
		}
		return doc;
	}

	public ContextDocument RemoveCitation(TIntIntHashMap citedDocuments) {
		ContextDocument doc = new ContextDocument();
		doc.word_counts = this.word_counts;
		doc.wordSequence = this.wordSequence;
		doc.topic_counts = this.topic_counts;
		doc.topicAssignments = this.topicAssignments;
		doc.docId = this.docId;
		doc.docName = this.docName;
		doc.docLength = this.docLength;
		doc.docLabel = this.docLabel;
		doc.wordLength = this.wordLength;
		doc.sentenseLength = this.sentenseLength;
		int[] keys = this.SentenseBoundaries.keys();
		for (int i = 0; i < keys.length; i++) {
			Vector<Integer> sentences = new Vector<Integer>();
			for (int j = 0; j < this.SentenseBoundaries.get(i).size(); j++) {
				sentences.add(this.SentenseBoundaries.get(i).get(j));
			}
			doc.SentenseBoundaries.put(keys[i], sentences);
		}
		// put only those citations which are in training set
		for (int i = 0; i < this.sentenseLength; i++) {
			if (this.contextObject.contains(i)) {
				Vector<Integer> citations = new Vector<Integer>();
				for (int j = 0; j < this.contextObject.get(i).size(); j++) {
					if (citedDocuments.contains(this.contextObject.get(i)
							.get(j))) {
						citations.add(this.contextObject.get(i).get(j));
						doc.citationSet
								.put(this.contextObject.get(i).get(j), 1);
					}
				}
				if (citations.size() > 0)
					doc.contextObject.put(i, citations);

			}
		}
		// copy the citationset in a seperate variable
		keys = this.citationSet.keys();
		for (int i = 0; i < keys.length; i++) {
			doc.citationSetTest.adjustOrPutValue(keys[i], 1, 1);
		}
		return doc;
	}

	public void setWindowLenghtAdaptive(Corpus corpus, int fold) {
		// for each citation context, get the top 10 lines of citation and its
		// topic count.
		int total_contexts = 0;
		for (int i = 0; i < citationContexts.size(); i++) {
			CitationContext cc = citationContexts.get(i);
			if (cc.length >= 10) {
				continue;
			}
			ContextDocument doc = (ContextDocument) corpus.docsTable
					.get(cc.citationid);
			if (doc==null || doc.topicAssignments == null) {
				continue;
			}
			total_contexts++;
			TIntDoubleHashMap citation_topics = new TIntDoubleHashMap();
			TIntDoubleHashMap context_topics = new TIntDoubleHashMap();
			TIntDoubleHashMap context_topics_left = new TIntDoubleHashMap();
			TIntDoubleHashMap context_topics_right = new TIntDoubleHashMap();
			TIntDoubleHashMap context_topics_both = new TIntDoubleHashMap();
			int[] topics = doc.topicAssignments.getFeatures();
			for (int l = 0; l < doc.sentenseLength && l < 10; l++) {
				for (int j = doc.SentenseBoundaries.get(l).firstElement(); j < doc.SentenseBoundaries
						.get(l).lastElement(); j++) {
					citation_topics.adjustOrPutValue(topics[j], 1, 1);
				}
			}
			topics = this.topicAssignments.getFeatures();
			for (int l = cc.begin; l <= cc.end; l++) {
				for (int j = this.SentenseBoundaries.get(l).firstElement(); j < this.SentenseBoundaries
						.get(l).lastElement(); j++) {
					context_topics.adjustOrPutValue(topics[j], 1, 1);

				}
			}
			int radius = 1;
			int left = cc.begin - radius, right = cc.end + radius;
			if (left < 0) {
				left = cc.begin;
			}
			if (right >= this.sentenseLength) {
				right = this.sentenseLength - 1;
			}
			for (int l = left; l <= cc.end; l++) {
				for (int j = this.SentenseBoundaries.get(l).firstElement(); j < this.SentenseBoundaries
						.get(l).lastElement(); j++) {

					context_topics_left.adjustOrPutValue(topics[j], 1, 1);

				}
			}
			// System.out.println(this.sentenseLength + ":"+right );
			for (int l = cc.begin; l <= right; l++) {
				for (int j = this.SentenseBoundaries.get(l).firstElement(); j < this.SentenseBoundaries
						.get(l).lastElement(); j++) {

					context_topics_right.adjustOrPutValue(topics[j], 1, 1);

				}
			}
			for (int l = left; l <= right; l++) {
				for (int j = this.SentenseBoundaries.get(l).firstElement(); j < this.SentenseBoundaries
						.get(l).lastElement(); j++) {

					context_topics_both.adjustOrPutValue(topics[j], 1, 1);
				}
			}
			citation_topics = normalizeSparseVector(citation_topics);
			context_topics = normalizeSparseVector(context_topics);
			context_topics_left = normalizeSparseVector(context_topics_left);
			context_topics_right = normalizeSparseVector(context_topics_right);
			context_topics_both = normalizeSparseVector(context_topics_both);
			// get similarity
			double max = topicVectorSimilarity(citation_topics, context_topics);
			int side = 0;
			double tmp = topicVectorSimilarity(citation_topics,
					context_topics_left);
			if (max < tmp) {
				side = 1;
				max = tmp;
			}
			tmp = topicVectorSimilarity(citation_topics, context_topics_right);
			if (max < tmp) {
				side = 2;
				max = tmp;
			}
			tmp = topicVectorSimilarity(citation_topics, context_topics_both);
			if (max < tmp) {
				side = 3;
				max = tmp;
			}

			// change context object for this citation and context.
			switch (side) {

			case 1:
				if (cc.begin > 0) {
					cc.begin--;
					avgContextLength -= cc.length;
					cc.updateLength();
					avgContextLength += cc.length;
					avgLeftContextLength++;

					if (!this.contextObject.contains(cc.begin)) {
						Vector<Integer> citations = new Vector<Integer>();
						citations.add(cc.citationid);
						this.contextObject.put(cc.begin, citations);

					} else {
						this.contextObject.get(cc.begin).add(cc.citationid);
					}
				}

				// increase the left side
				break;
			case 2:
				if (cc.end < this.sentenseLength - 1) {
					cc.end++;
					avgContextLength -= cc.length;
					cc.updateLength();
					avgContextLength += cc.length;
					avgRightContextLength++;
					if (!this.contextObject.contains(cc.end)) {
						Vector<Integer> citations = new Vector<Integer>();
						citations.add(cc.citationid);
						this.contextObject.put(cc.end, citations);

					} else {
						this.contextObject.get(cc.end).add(cc.citationid);
					}
				}

				// increase the right side
				break;
			case 3:
				if (cc.begin > 0) {
					cc.begin--;
					avgContextLength -= cc.length;
					cc.updateLength();
					avgContextLength += cc.length;
					avgLeftContextLength++;
					if (!this.contextObject.contains(cc.begin)) {
						Vector<Integer> citations = new Vector<Integer>();
						citations.add(cc.citationid);
						this.contextObject.put(cc.begin, citations);

					} else {
						this.contextObject.get(cc.begin).add(cc.citationid);
					}
				}
				if (cc.end < this.sentenseLength - 1) {
					cc.end++;
					avgContextLength -= cc.length;
					cc.updateLength();
					avgContextLength += cc.length;
					avgRightContextLength++;
					if (!this.contextObject.contains(cc.end)) {
						Vector<Integer> citations = new Vector<Integer>();
						citations.add(cc.citationid);
						this.contextObject.put(cc.end, citations);

					} else {
						this.contextObject.get(cc.end).add(cc.citationid);
					}
				}

			default:
				break;
			}
		}
		if (total_contexts > 0) {
			avgContextLength /= total_contexts;
			avgLeftContextLength /= total_contexts;
			avgRightContextLength /= total_contexts;
		}

	}

}