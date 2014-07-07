package edu.psu.types;

import java.util.*;

import edu.psu.types.*;
import edu.psu.util.ReadDirectory;
import gnu.trove.*;

import java.io.*;
import java.math.*;

public class Corpus implements Serializable {
	public static Alphabet vocabulary = new Alphabet();
	public static Alphabet citationAlphabet = new Alphabet();
	public static Alphabet docAlphabet = new Alphabet();
	public TIntDoubleHashMap[] cluster_states;
	public TIntDoubleHashMap[] cluster_probs_reconst;
	public static double[] mean;
	public double[][] eignVector;
	public double[][] normalizedEignVector;
	public double[][][] subspaceEignVector;
	public double[] mean_projection;
	public double[][] c_proj;
	public double[][] cluster;
	public int[] cluster_doc;
	public int[] cluster_size;
	public double[][] belong;
	public double[] cluster_probs;
	public TIntIntHashMap[] word_class_occ;

	public double[] alphas;
	public TIntDoubleHashMap[] subspace_vectors;
	public int[] doc_cluster;
	public int vocab_size;
	public static int size;
	public Vector<Document> docs;
	public TIntObjectHashMap<Document> docsTable;
	public Vector<WordOccurrences> words;
	public static TIntDoubleHashMap idf;
	public TIntIntHashMap[] typeTopicCounts;
	public TIntIntHashMap[] citationTopicCounts;
	public double beta;
	public double betaSum;
	public double gamma;
	public double gammaSum;
	public double[] alpha;
	public double alphaSum;
	public int numTypes;
	public int numCitations;
	public int[] tokensPerTopic;
	public int[] citationsPerTopic;
	public double[][] theta_train;
	public double[][] phi_train;
	public double[][] psi_train;
	public boolean isPerplex = false;
	public static TIntIntHashMap split;
	public static Map<Integer, Vector<Document>> test_docs = new HashMap<Integer, Vector<Document>>();
	public static Map<Integer, Vector<Document>> train_docs = new HashMap<Integer, Vector<Document>>();
	public static Map<Integer, ConceptVectors> concepts = new HashMap<Integer, ConceptVectors>();
	public TIntIntHashMap[] trainCitations;
	public Vector<RankedList> citedTopicLists;
	public static RankedList citedList = new RankedList();
	public int maxTokens;

	public Corpus() {

	}

	public class Splits {
		TIntIntHashMap[] tests;
		TIntIntHashMap[] trains;

		public Splits(int numFolds) {
			tests = new TIntIntHashMap[numFolds];
			trains = new TIntIntHashMap[numFolds];
			for (int i = 0; i < numFolds; i++) {
				tests[i] = new TIntIntHashMap();
				trains[i] = new TIntIntHashMap();
			}
		}

		public void produceSplits(Integer[] instances, int numFolds) {
			// shuffle array
			Random r = new Random(21);
			System.out.println(instances.length);
			for (int i = 0; i < instances.length; i++) {
				int rand = r.nextInt(instances.length - i) + i;
				int temp = instances[i];
				instances[i] = instances[rand];
				instances[rand] = temp;
			}
			for (int i = 0; i < instances.length; i++) {
				int which_split = i % numFolds;
				// System.out.println(instances[i].intValue());

				tests[which_split].put(instances[i].intValue(), 1);
				for (int j = 0; j < numFolds; j++) {
					if (j != which_split) {
						trains[j].put(instances[i].intValue(), 1);
					}
				}

			}

		}

	}

	public void prepareSplits(int numFolds, boolean isCitationInformation) {

		Integer[] entries = new Integer[docs.size()];
		for (int i = 0; i < docs.size(); i++) {
			entries[i] = i;
		}
		trainCitations = new TIntIntHashMap[numFolds];
		for (int i = 0; i < numFolds; i++) {
			Vector<Document> test = new Vector<Document>();
			Vector<Document> train = new Vector<Document>();
			test_docs.put(i, test);
			train_docs.put(i, train);
		}
		Splits split = new Splits(numFolds);
		split.produceSplits(entries, numFolds);
		for (int i = 0; i < numFolds; i++) {
			trainCitations[i] = new TIntIntHashMap();
			int[] keys = split.trains[i].keys();
			for (int j = 0; j < keys.length; j++) {
				ContextDocument doc = (ContextDocument) this.docs.get(keys[j]);
				int[] cited = doc.citationSet.keys();
				for (int k = 0; k < cited.length; k++) {
					trainCitations[i].put(cited[k], 1);
				}

			}
			keys = split.tests[i].keys();
			System.out
					.println("Before removing doc with no cited articles in fold"
							+ i
							+ ": train="
							+ split.trains[i].size()
							+ ":Test=" + split.tests[i].size());
			// take documents away from test set which dont have citations.
			for (int j = 0; j < keys.length; j++) {
				ContextDocument doc = (ContextDocument) this.docs.get(keys[j]);
				if (doc.citationSet.size() == 0) {
					split.trains[i].put(keys[j], 1);
					split.tests[i].remove(keys[j]);
				} else {
					int[] cited = doc.citationSet.keys();
					boolean seen_in_training = false;
					for (int k = 0; k < cited.length; k++) {
						if (trainCitations[i].containsKey(cited[k])) {
							seen_in_training = true;
						}
					}
					if (!seen_in_training) {
						split.trains[i].put(keys[j], 1);
						split.tests[i].remove(keys[j]);
					}
				}
			}
			System.out
					.println("After removing doc with no cited articles in fold"
							+ i
							+ ": train="
							+ split.trains[i].size()
							+ ":Test=" + split.tests[i].size());

		}

		// remove all the citation from test which do not appear in training.

		for (int i = 0; i < numFolds; i++) {
			Vector<Document> test = test_docs.get(i);
			Vector<Document> train = train_docs.get(i);
			int[] Keys = split.tests[i].keys();
			for (int j = 0; j < Keys.length; j++) {
				ContextDocument doc = (ContextDocument) docs.get(Keys[j]);
				test.add(doc.RemoveCitation(trainCitations[i]));
			}
			Keys = split.trains[i].keys();
			for (int j = 0; j < Keys.length; j++) {
				train.add(docs.get(Keys[j]));
			}
			// remove citation information
			// clone the documents in test and train splits and
			// remove citation information from test documents
			// do for citation recommendation only--check to see if test
			// document has any outgoing citation link or not

		}
		if (!isCitationInformation) {
			for (int i = 0; i < numFolds; i++) {
				Vector<Document> test = test_docs.get(i);

				for (int j = 0; j < test.size(); j++) {
					// System.out.print(((ContextDocument)test.get(j)).docId+":"+((ContextDocument)test.get(j)).citationSet.size()+":");
					ContextDocument doc = (ContextDocument) test.get(j);

					doc = doc.RemoveCitation();

					test.remove(j);
					test.insertElementAt(doc, j);
					// System.out.println(((ContextDocument)test.get(j)).docId+":"+((ContextDocument)test.get(j)).citationSet.size()+":");
				}
				// remove citation information
				// clone the documents in test and train splits and
				// remove citation information from test documents
				// do for citation recommendation only--check to see if test
				// document has any outgoing citation link or not

			}
		}

	}

	public static Corpus assignDocumentIds(String directory) {
		String[] files = ReadDirectory.list(directory);
		System.out.println("There are total " + files.length + " Files in "
				+ directory);
		for (int i = 0; i < files.length; i++) {
			String name = files[i].split("/")[files[i].split("/").length - 1]
					.replace(".txt", "");// check!!!
			if (!Corpus.docAlphabet.contains(name)) {
				Corpus.docAlphabet.lookupIndex(name);
			}
		}
		Corpus corpus = new Corpus(Corpus.docAlphabet.size());
		return corpus;

	}

	public void readData(String directory) {
		String[] files = ReadDirectory.list(directory);
		System.out.println("There are total " + files.length + " Files in "
				+ directory);
		String line = null;
		maxTokens=0;

		try {
			for (int i = 0; i < files.length; i++) {
				String name = files[i].split("/")[files[i].split("/").length - 1]
						.replace(".txt", "");// check!!!
				FileReader fr = new FileReader(files[i]);
				BufferedReader br = new BufferedReader(fr);
				try {
					Vector<String> lines = new Vector<String>();
					while ((line = br.readLine()) != null) {
						lines.add(line);
					}
					ContextDocument doc = new ContextDocument();
					doc.add(name, Corpus.docAlphabet.lookupIndex(name), lines);
					this.addDocument(doc, Corpus.docAlphabet.lookupIndex(name));
					if (maxTokens < doc.wordLength)
						maxTokens = doc.wordLength;
					fr.close();
					br.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

			}
			for (int i = 0; i < docsTable.size(); i++) {
				ContextDocument doc = (ContextDocument) docsTable.get(i);
				for (int j = 0; j < doc.citationSet.size(); j++) {
					Corpus.citedList.add(doc.citationSet.get(j), 1);
				}
			}

			FileWriter fw = new FileWriter("./cited_docs");
			BufferedWriter bw = new BufferedWriter(fw);

			for (int i = 0; i < Corpus.citedList.indices.size(); i++) {
				bw.write(Corpus.docAlphabet
						.lookupObject(Corpus.citedList.indices.get(i))
						+ ":"
						+ Corpus.citedList.values.get(Corpus.citedList.indices
								.get(i)) + "\n");
			}
			bw.close();
			fw.close();

			System.out.println("Total Cited Documents:"
					+ Corpus.citationAlphabet.size());
			System.out.println("Total Documents:" + this.docs.size());
			System.out.println("Total Vcabulary Size:"
					+ Corpus.vocabulary.size());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void setWindowLength(int length) {
		for (int i = 0; i < this.docs.size(); i++) {
			ContextDocument doc = (ContextDocument) this.docs.get(i);

			TIntObjectHashMap<Vector<Integer>> contextObject = doc.contextObject;
			int[] keys = contextObject.keys();
			for (int k = 0; k < keys.length; k++) {
				Vector<Integer> tmpContext = contextObject.get(keys[k]);
				for (int j = 1; j <= length; j++) {
					if (keys[k] + j < doc.sentenseLength) {
						contextObject.put(keys[k] + j,
								(Vector<Integer>) tmpContext.clone());
					}
					if (keys[k] - j >= 0) {
						contextObject.put(keys[k] - j,
								(Vector<Integer>) tmpContext.clone());
					}
				}
			}
		}
	}

	public void adaptiveWindowLength(int fold) {
		// assimilate topic assignment from present context, context+1,
		// context-1 and cited document's initial 10 sentences.
		System.out.print("Updating the context windows.......");
		int total = 0;
		double avgContextLength = 0, avgLeftContextLength = 0, avgRightContextLength = 0;
		for (int i = 0; i < this.docs.size(); i++) {
			ContextDocument doc = (ContextDocument) this.docs.get(i);
			doc.setWindowLenghtAdaptive(this, fold);
			if (doc.avgContextLength > 0) {
				total++;
			}
			avgContextLength += doc.avgContextLength;
			avgLeftContextLength += doc.avgLeftContextLength;
			avgRightContextLength += doc.avgRightContextLength;
		}
		avgContextLength /= total;
		avgLeftContextLength /= total;
		avgRightContextLength /= total;
		System.out.println("Done!:Avg Context Length=" + avgContextLength
				+ ":Avg Left Context Length=" + avgLeftContextLength
				+ ":Avg Right Context Length=" + avgRightContextLength);
	}

	public Corpus(int capacity) {
		docs = new Vector<Document>();
		for (int i = 0; i < capacity; i++) {
			docs.add(new Document());
		}
		size = capacity;
	}

	public void addDocument(Document doc) {
		if (docs == null) {
			docs = new Vector<Document>();
		}
		docs.add(doc);
		size = docs.size();
	}

	public void addDocument(Document doc, int index) {
		docs.set(index, doc);
		if (docsTable == null) {
			docsTable = new TIntObjectHashMap<Document>();
		}
		docsTable.put(index, doc);
	}

	public Document getDocument(int index) {
		if (docs.size() <= index) {
			System.err.println("doc id exceed corpus size!");
			System.exit(-1);
		}
		return docs.get(index);
	}

	public void addDocuments(Vector<Document> documents) {

		docs = documents;
		size = docs.size();
	}

	public void prepareWordOccurrences() {
		words = new Vector<WordOccurrences>();
		for (int i = 0; i < this.vocabulary.size(); i++) {
			words.add(new WordOccurrences(i, this));
		}
	}

	public void idf() {
		if (docs == null || docs.size() == 0) {
			System.out.println("No Documents. So can not prepare idf vector");
			idf = null;
			return;
		}
		idf = new TIntDoubleHashMap();
		TIntIntHashMap tmp = new TIntIntHashMap();

		// prepare word class occurrences
		word_class_occ = new TIntIntHashMap[2];
		for (int i = 0; i < 2; i++)
			word_class_occ[i] = new TIntIntHashMap();
		for (Document doc : docs) {
			String name = doc.docName;
			String labl = name.split("/")[0];

			int[] wordids = doc.word_counts.keys();
			for (int i = 0; i < wordids.length; i++) {
				if (labl.contains("guns")) {
					word_class_occ[0].adjustOrPutValue(wordids[i], 1, 1);
				} else
					word_class_occ[1].adjustOrPutValue(wordids[i], 1, 1);
				tmp.adjustOrPutValue(wordids[i], 1, 1);
			}
		}
		int vocab_size = vocabulary.size();
		int corpus_size = docs.size();
		for (int i = 0; i < vocab_size; i++) {
			double val = Math.log((double) corpus_size
					/ (double) (tmp.get(i) + 1));
			idf.adjustOrPutValue(i, val, val);
		}

	}

	public void normalizeCorpus() {
		this.idf();
		this.setMean();
		for (Document doc : docs) {
			doc.normalizeTFIDF();
		}
	}

	public void setMean() {
		mean = new double[this.vocabulary.size()];
		for (int i = 0; i < docs.size(); i++) {
			TIntIntHashMap doc = docs.get(i).word_counts;
			int[] keys = doc.keys();
			for (int j = 0; j < keys.length; j++) {
				mean[keys[j]] += (doc.get(keys[j]) / (double) docs.size());
			}
		}

		/*
		 * double norm = 0.0; for (int i = 0; i < mean.length; i++) { norm +=
		 * Math.pow(mean[i], 2); } norm = Math.pow(norm, 0.5); for (int i = 0; i
		 * < mean.length; i++) { mean[i] = mean[i] / norm; }
		 */
		System.out.println("length=" + mean.length);

	}

	public void normalizeClusterCentroid(int cluster, double p) {

		double norm = 0.0;
		int[] words = cluster_states[cluster].keys();
		for (int j = 0; j < words.length; j++) {
			norm += Math.pow(cluster_states[cluster].get(words[j]), p);
		}
		norm = Math.pow(norm, 1.0 / p);
		for (int j = 0; j < words.length; j++) {
			if (norm == 0) {
				System.out.println(j);
			} else
				cluster_states[cluster].put(words[j],
						cluster_states[cluster].get(words[j]) / norm);
			// System.out.println(cluster_states[i].get(words[j]));
		}

	}

	public void normalizeSubspaceCentroid(int cluster, double p) {

		double norm = 0.0;
		int[] words = this.subspace_vectors[cluster].keys();
		for (int j = 0; j < words.length; j++) {
			norm += Math.pow(this.subspace_vectors[cluster].get(words[j]), p);
		}
		norm = Math.pow(norm, 1.0 / p);
		for (int j = 0; j < words.length; j++) {
			if (norm == 0) {
				System.out.println(j);
			} else
				this.subspace_vectors[cluster].put(words[j],
						this.subspace_vectors[cluster].get(words[j]) / norm);
			// System.out.println(cluster_states[i].get(words[j]));
		}

	}

	public void normalizeClusterCentroids(double p) {
		for (int i = 0; i < cluster_states.length; i++) {
			double norm = 0.0;
			int[] words = cluster_states[i].keys();
			for (int j = 0; j < words.length; j++) {
				norm += Math.pow(cluster_states[i].get(words[j]), p);
			}
			norm = Math.pow(norm, 1.0 / p);
			for (int j = 0; j < words.length; j++) {
				if (norm == 0) {
					System.out.println(j);
				} else
					cluster_states[i].put(words[j],
							cluster_states[i].get(words[j]) / norm);
				// System.out.println(cluster_states[i].get(words[j]));
			}

		}
	}

	public void printTopWordsSubSpace(BufferedWriter output, int num_clusters,
			int count, boolean reconstruct) {
		try {
			for (int i = 0; i < num_clusters; i++) {

				TreeSet<IDSorter> sortedWords = new TreeSet<IDSorter>();
				int[] keys = this.subspace_vectors[i].keys();
				for (int k = 0; k < keys.length; k++) {
					sortedWords.add(new IDSorter(keys[k],
							this.subspace_vectors[i].get(keys[k])));
				}
				output.write("Cluster" + i + ":prob=" + this.alphas[i]);
				output.write("\n------------\n");
				Iterator<IDSorter> iterator = sortedWords.iterator();
				int words = 1;
				// while(iterator.hasNext() && words<keys.length-count){
				// IDSorter info = iterator.next();
				// words++;
				// }

				while (iterator.hasNext() && words < count) {
					IDSorter info = iterator.next();
					// output.write(this.vocabulary.lookupObject(info.getID())
					// + ":" + info.getWeight() + ":"
					// + this.cluster_states[i].get(info.getID()) + "\n");
					output.write(this.vocabulary.lookupObject(info.getID())
							+ ":" + info.getWeight() + "\t");
					words++;
				}
				// System.out.println("-----"+words+"-------");
				output.write("\n------------\n");

			}
			output.write("\nCluster Probabilities:\n\n");
			for (int i = 0; i < num_clusters; i++) {
				output.write(i + "=" + this.alphas[i] + "\n");
			}
			output.write("------------\n\n");
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

	}

	public void printTopWordsEign(BufferedWriter output, int num_clusters,
			int count) {
		try {
			double[][] cluster1 = new double[5][this.vocabulary.size()];

			TreeSet<IDSorter>[] SortedWords = new TreeSet[this.normalizedEignVector.length];
			for (int i = 0; i < SortedWords.length; i++) {
				SortedWords[i] = new TreeSet<IDSorter>();
				for (int k = 0; k < vocabulary.size(); k++) {
					SortedWords[i].add(new IDSorter(k,
							this.normalizedEignVector[i][k]));
				}
			}
			Iterator<IDSorter>[] iterator = new Iterator[this.normalizedEignVector.length];
			for (int i = 0; i < SortedWords.length; i++) {
				iterator[i] = SortedWords[i].iterator();
			}
			int words = 0;

			while (iterator[0].hasNext()) {
				IDSorter info = iterator[0].next();
				if (words < count || words > this.vocabulary.size() - count) {
					output.write(this.vocabulary.lookupObject(info.getID())
							+ ":" + info.getWeight() + "\t");
				}

				for (int i = 1; i < SortedWords.length; i++) {
					info = iterator[i].next();
					if (words < count || words > this.vocabulary.size() - count) {
						output.write(this.vocabulary.lookupObject(info.getID())
								+ ":" + info.getWeight() + "\t");
					}
				}
				if (words <= count || words >= this.vocabulary.size() - count) {
					if (words == this.vocabulary.size() - count)
						output.write("------------------------------------------------------------\n");
					else
						output.write("\n");
				}
				words++;

			}

			output.write("\n------------\n");

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

	}

	public void printTopWords(BufferedWriter output, int num_clusters,
			int count, boolean reconstruct) {
		if (!reconstruct) {
			try {
				for (int i = 0; i < num_clusters; i++) {

					TreeSet<IDSorter> sortedWords = new TreeSet<IDSorter>();
					int[] keys = this.cluster_states[i].keys();
					for (int k = 0; k < keys.length; k++) {
						sortedWords.add(new IDSorter(keys[k],
								this.cluster_states[i].get(keys[k])));
					}
					output.write("Cluster" + i + ":" + keys.length + ":prob="
							+ this.cluster_probs[i]);
					output.write("\n------------\n");
					Iterator<IDSorter> iterator = sortedWords.iterator();
					int words = 1;
					// while(iterator.hasNext() && words<keys.length-count){
					// IDSorter info = iterator.next();
					// words++;
					// }

					while (iterator.hasNext() && words < count) {
						IDSorter info = iterator.next();
						// output.write(this.vocabulary.lookupObject(info.getID())
						// + ":" + info.getWeight() + ":"
						// + this.cluster_states[i].get(info.getID()) + "\n");
						output.write(this.vocabulary.lookupObject(info.getID())
								+ ":" + info.getWeight() + "\t");
						words++;
					}
					// System.out.println("-----"+words+"-------");
					output.write("\n------------\n");

				}
				output.write("\nCluster Probabilities:\n\n");
				for (int i = 0; i < num_clusters; i++) {
					output.write(i + "=" + this.cluster_probs[i] + "\n");
				}
				output.write("------------\n\n");
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(0);
			}
		} else {
			try {
				for (int i = 0; i < num_clusters; i++) {

					TreeSet<IDSorter> sortedWords = new TreeSet<IDSorter>();
					int[] keys = this.cluster_states[i].keys();
					for (int k = 0; k < keys.length; k++) {
						sortedWords.add(new IDSorter(keys[k],
								this.cluster_states[i].get(keys[k])));
					}
					output.write("Cluster" + i + ":");

					Iterator<IDSorter> iterator = sortedWords.iterator();
					int words = 1;
					// while(iterator.hasNext() && words<keys.length-count){
					// IDSorter info = iterator.next();
					// words++;
					// }

					while (iterator.hasNext() && words < count) {
						IDSorter info = iterator.next();
						// output.write(this.vocabulary.lookupObject(info.getID())
						// + ":" + info.getWeight() + ":"
						// + this.cluster_states[i].get(info.getID()) + "\n");
						output.write(this.vocabulary.lookupObject(info.getID())
								+ ":" + info.getWeight() + "\t");
						words++;
					}
					// System.out.println("-----"+words+"-------");
					output.write("\n------------\n");

				}

				output.write("------------\n\n");
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(0);
			}
		}
	}

	public void printDocBelongings(BufferedWriter output, int num_clusters) {
		try {
			for (int i = 0; i < this.docs.size(); i++) {
				TreeSet<IDSorter> sortedWords = new TreeSet<IDSorter>();
				for (int j = 0; j < num_clusters; j++) {
					sortedWords.add(new IDSorter(j,
							this.cluster_probs_reconst[j].get(i)));
				}
				Iterator<IDSorter> iterator = sortedWords.iterator();
				output.write(this.docs.get(i).docName + " ::: ");
				while (iterator.hasNext()) {
					IDSorter info = iterator.next();
					output.write(info.getID() + ":" + info.getWeight() + "\t");
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	public void printFrequentWords(BufferedWriter output,
			TIntDoubleHashMap vector, int count) {

		try {
			TreeSet<IDSorter> sortedWords = new TreeSet<IDSorter>();
			int[] keys = vector.keys();
			for (int k = 0; k < keys.length; k++) {
				sortedWords.add(new IDSorter(keys[k], vector.get(keys[k])));
			}

			Iterator<IDSorter> iterator = sortedWords.iterator();
			int words = 1;
			// while(iterator.hasNext() && words<keys.length-count){
			// IDSorter info = iterator.next();
			// words++;
			// }

			while (iterator.hasNext() && words < count) {
				IDSorter info = iterator.next();
				output.write(this.vocabulary.lookupObject(info.getID()) + ":"
						+ info.getWeight() + "\n");
				words++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		// System.out.println("-----"+words+"-------");

	}

	public void printFrequentWords(BufferedWriter output, double[] vector,
			int count) {

		try {
			TreeSet<IDSorter> sortedWords = new TreeSet<IDSorter>();

			for (int k = 0; k < vector.length; k++) {
				sortedWords.add(new IDSorter(k, vector[k]));
			}

			Iterator<IDSorter> iterator = sortedWords.iterator();
			int words = 1;
			// while(iterator.hasNext() && words<keys.length-count){
			// IDSorter info = iterator.next();
			// words++;
			// }

			while (iterator.hasNext() && words < count) {
				IDSorter info = iterator.next();
				output.write(this.vocabulary.lookupObject(info.getID()) + ":"
						+ info.getWeight() + "\n");
				words++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		// System.out.println("-----"+words+"-------");

	}

	public void printFrequentWords(TIntDoubleHashMap vector, int count) {

		TreeSet<IDSorter> sortedWords = new TreeSet<IDSorter>();
		int[] keys = vector.keys();
		for (int k = 0; k < keys.length; k++) {
			sortedWords.add(new IDSorter(keys[k], vector.get(keys[k])));
		}

		Iterator<IDSorter> iterator = sortedWords.iterator();
		int words = 1;
		// while(iterator.hasNext() && words<keys.length-count){
		// IDSorter info = iterator.next();
		// words++;
		// }

		while (iterator.hasNext() && words < count) {
			IDSorter info = iterator.next();
			System.out.println(this.vocabulary.lookupObject(info.getID()) + ":"
					+ info.getWeight());
			words++;
		}
		// System.out.println("-----"+words+"-------");
		System.out.println("------------");

	}

	public void printPairwiseDist(BufferedWriter output) {

		try {
			output.write("Pairwise Cluster Distances:\n\n");
			double dotp;
			for (int j = 0; j < cluster_states.length; j++) {
				for (int i = 0; i < cluster_states.length; i++) {

					dotp = 0.0;
					int[] words = cluster_states[i].keys();

					for (int k = 0; k < words.length; k++) {
						if (this.cluster_states[j].containsKey(words[k]))
							dotp += this.cluster_states[j].get(words[k])
									* this.cluster_states[i].get(words[k]);
					}
					output.write(j + ":" + i + "=" + dotp + "\n");

				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	public static void normalizeVector(double[] vector, double p) {
		double norm = 0.0;
		for (int i = 0; i < vector.length; i++) {
			norm += Math.pow(vector[i], p);
		}
		norm = Math.pow(norm, 1.0 / p);
		for (int i = 0; i < vector.length; i++) {
			vector[i] /= norm;
		}
	}

	public double recostruction_error(int num_clusters) {
		double err = 0;
		int vocab = this.vocabulary.size();
		for (int i = 0; i < this.docs.size(); i++) {
			for (int v1 = 0; v1 < vocab; v1++) {
				for (int v2 = 0; v2 < vocab; v2++) {
					// err+=Math.pow((this.c_proj[v1][v2]*-),2);

				}

			}

		}

		return err;

	}

	public double loglikelihood(int num_clusters) {
		double ll = 0.0, dotp, cum_resp;
		// for(int j=0;j<num_clusters;j++)
		// System.out.println(this.cluster_probs[j]);
		for (int i = 0; i < this.docs.size(); i++) {

			cum_resp = 0;
			int[] keys_data = this.docs.get(i).normed_tfidfVector.keys();
			for (int j = 0; j < num_clusters; j++) {
				dotp = 0;
				for (int k = 0; k < keys_data.length; k++) {
					if (this.cluster_states[j].containsKey(keys_data[k]))
						dotp += this.docs.get(i).normed_tfidfVector
								.get(keys_data[k])
								* this.cluster_states[j].get(keys_data[k]);
				}
				cum_resp += dotp * dotp * this.cluster_probs[j];

			}
			if (cum_resp > 0)
				ll += Math.log(cum_resp);
		}
		return ll;

	}

	public void printSubspaceDocRelevance(BufferedWriter output,
			int num_clusters, int count) {
		int vocab = this.vocabulary.size();
		try {
			for (int i = 0; i < num_clusters; i++) {
				TreeSet<IDSorter> sortedWords = new TreeSet<IDSorter>();
				for (int j = 0; j < this.docs.size(); j++) {
					double[] inter_vec = new double[vocab];
					double prob = 0;
					TIntDoubleHashMap doc = this.docs.get(j).normed_tfidfVector;
					int[] keys = doc.keys();
					for (int v = 0; v < vocab; v++) {
						for (int k = 0; k < keys.length; k++) {
							inter_vec[v] += doc.get(keys[k])
									* this.c_proj[keys[k]][v];
						}
					}
					for (int k = 0; k < keys.length; k++) {
						prob += doc.get(keys[k]) * inter_vec[keys[k]];
					}
					sortedWords.add(new IDSorter(j, prob));
				}
				TObjectIntHashMap<String> lablCount = new TObjectIntHashMap<String>();
				Iterator<IDSorter> iterator = sortedWords.iterator();

				int docs = 1;
				// while(iterator.hasNext() && words<keys.length-count){
				// IDSorter info = iterator.next();
				// words++;
				// }
				output.write("\n------------\n");
				output.write("Cluster-" + i + ": Prob:" + this.cluster_probs[i]
						+ "\n");
				while (iterator.hasNext() && docs < 500) {
					IDSorter info = iterator.next();
					String name = this.docs.get(info.getID()).docName;
					String labl = name.split("/")[0];
					lablCount.adjustOrPutValue(labl, 1, 1);
					// if (docs < count)
					// output.write(this.docs.get(info.getID()).docName + ":"
					// + info.getWeight() + "\t");
					docs++;

				}
				output.write("\n");
				Object[] keys = lablCount.keys();

				for (Object key : keys) {
					output.write((String) key + ":"
							+ lablCount.get((String) key) + "\t");
				}

			}

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	/*
	 * [a b][i j][a] = a^2*i+a*b*k+a*b*j+b*2*l [k l][b]
	 */
	public void printRankSubDocRelevance(BufferedWriter output,
			int num_clusters, int count) {
		try {
			for (int i = 0; i < num_clusters; i++) {
				TreeSet<IDSorter> sortedWords = new TreeSet<IDSorter>();
				for (int j = 0; j < this.docs.size(); j++) {
					int[] keys = this.docs.get(j).normed_tfidfVector.keys();
					double dotp = 0;
					for (int k = 0; k < keys.length; k++) {
						dotp += this.docs.get(j).normed_tfidfVector
								.get(keys[k]) * this.c_proj[i][keys[k]];
						output.write(this.c_proj[i][keys[k]] + "\n");
					}

					sortedWords.add(new IDSorter(j, dotp));

				}
				TObjectIntHashMap<String> lablCount = new TObjectIntHashMap<String>();
				Iterator<IDSorter> iterator = sortedWords.iterator();

				int docs = 1;
				// while(iterator.hasNext() && words<keys.length-count){
				// IDSorter info = iterator.next();
				// words++;
				// }
				// output.write("\n------------\n");
				// output.write("Cluster-" + i + ": Prob:" +
				// this.cluster_probs[i]
				// + "\n");
				while (iterator.hasNext() && docs < 500) {
					IDSorter info = iterator.next();
					String name = this.docs.get(info.getID()).docName;
					String labl = name.split("/")[0];
					lablCount.adjustOrPutValue(labl, 1, 1);
					if (docs < count)
						output.write(this.docs.get(info.getID()).docName + ":"
								+ info.getWeight() + "\t");
					docs++;

				}
				output.write("\n");
				Object[] keys = lablCount.keys();

				for (Object key : keys) {
					output.write((String) key + ":"
							+ lablCount.get((String) key) + "\t");
				}
				// System.out.println("-----"+words+"-------");

			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

	}

	public void printClustersNMF(BufferedWriter output, int num_clusters) {
		try {
			int[] cluster_size = new int[num_clusters];
			TObjectIntHashMap<String>[] class_clusters = new TObjectIntHashMap[num_clusters];
			for (int i = 0; i < num_clusters; i++) {
				class_clusters[i] = new TObjectIntHashMap<String>();
			}

			for (int j = 0; j < this.docs.size(); j++) {

				String name = this.docs.get(j).docName;
				String labl = name.split("/")[0];
				class_clusters[this.cluster_doc[j]]
						.adjustOrPutValue(labl, 1, 1);
				cluster_size[this.cluster_doc[j]]++;

			}
			output.write("\n\n-----------Document Distribution in clusters--------------\n");
			for (int i = 0; i < num_clusters; i++) {

				TreeSet<DescSorter> sortedCat = new TreeSet<DescSorter>();
				Object[] keys = class_clusters[i].keys();
				for (Object key : keys) {
					sortedCat.add(new DescSorter((String) key,
							class_clusters[i].get((String) key)));
				}
				Iterator<DescSorter> DescIter = sortedCat.iterator();
				output.write("Cluster-" + i + ": # docs in cluster= "
						+ cluster_size[i] + ":\n");
				while (DescIter.hasNext()) {
					DescSorter info = DescIter.next();
					output.write(info.getID() + ":" + info.getWeight() + "\t");
				}
				output.write("\n\n");
			}
			output.write("\n\n-----------Cluster Entropy for sub class--------------\n");
			double entropy = 0;
			for (int i = 0; i < num_clusters; i++) {
				double cluster_entropy = 0;
				if (cluster_size[i] > 0) {
					for (Object key : class_clusters[i].keys()) {
						cluster_entropy += (class_clusters[i].get((String) key) * Math
								.log((double) class_clusters[i]
										.get((String) key)
										/ (double) cluster_size[i]));
						// System.out.println(class_clusters[i].get((String)
						// key)+":"+cluster_size[i]);
					}

					cluster_entropy = (-1) * cluster_entropy
							/ (double) this.docs.size();

				}
				output.write("Cluster-" + i + ": # docs in cluster= "
						+ cluster_size[i] + ":");

				output.write("Class Entropy :" + cluster_entropy);

				output.write("\n\n");
				entropy += cluster_entropy;
			}
			output.write("Overall Entropy :" + entropy + "\n\n");
			System.out.println("Overall Entropy :" + entropy + "\n\n");

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

	}

	public void printClusters(BufferedWriter output, int num_clusters, int count) {
		try {
			int[] cluster_size = new int[num_clusters];
			TObjectIntHashMap<String>[] class_clusters = new TObjectIntHashMap[num_clusters];
			for (int i = 0; i < num_clusters; i++) {
				class_clusters[i] = new TObjectIntHashMap<String>();
			}

			for (int j = 0; j < this.docs.size(); j++) {
				int max_cluster = -1;
				double max_dotp = 0;
				for (int i = 0; i < num_clusters; i++) {
					int[] keys = this.docs.get(j).normed_tfidfVector.keys();
					double dotp = 0;
					for (int k = 0; k < keys.length; k++) {
						dotp += this.docs.get(j).normed_tfidfVector
								.get(keys[k])
								* this.cluster_states[i].get(keys[k]);
					}
					if (dotp > max_dotp) {
						max_cluster = i;
						max_dotp = dotp;
					}

				}
				if (max_cluster != -1) {
					String name = this.docs.get(j).docName;
					String labl = name.split("/")[0];
					class_clusters[max_cluster].adjustOrPutValue(labl, 1, 1);
					cluster_size[max_cluster]++;
				} else {
					// System.out.println("Error in finding relevance cluster!!");
				}

			}
			output.write("\n\n-----------Document Distribution in clusters--------------\n");
			for (int i = 0; i < num_clusters; i++) {

				TreeSet<DescSorter> sortedCat = new TreeSet<DescSorter>();
				Object[] keys = class_clusters[i].keys();
				for (Object key : keys) {
					sortedCat.add(new DescSorter((String) key,
							class_clusters[i].get((String) key)));
				}
				Iterator<DescSorter> DescIter = sortedCat.iterator();
				output.write("Cluster-" + i + ": Prob:" + this.cluster_probs[i]
						+ ": # docs in cluster= " + cluster_size[i] + ":\n");
				while (DescIter.hasNext()) {
					DescSorter info = DescIter.next();
					output.write(info.getID() + ":" + info.getWeight() + "\t");
				}
				output.write("\n\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

	}

	public void printEignClusters(BufferedWriter output, int num_clusters,
			int count) {
		try {
			double[][] cluster1 = new double[5][this.vocabulary.size()];

			for (int i = 0; i < this.vocabulary.size(); i++) {
				cluster1[0][i] = (this.normalizedEignVector[0][i]);
				cluster1[1][i] = (this.normalizedEignVector[1][i]);
				cluster1[2][i] = (this.normalizedEignVector[2][i]);
				cluster1[3][i] = (this.normalizedEignVector[3][i]);
				cluster1[4][i] = (this.normalizedEignVector[4][i]);

			}
			for (int t = 0; t < 5; t++) {
				double norm = 0.0;
				for (int i = 0; i < cluster1[t].length; i++) {
					norm += Math.pow(cluster1[t][i], 2);
				}
				norm = Math.pow(norm, 0.5);
				for (int i = 0; i < mean.length; i++) {
					cluster1[0][i] = cluster1[t][i] / norm;
				}
			}
			this.mean_projection = new double[5];
			double dotp_clus = 0;
			for (int i = 0; i < this.vocabulary.size(); i++) {
				dotp_clus += cluster1[0][i] * cluster1[1][i];
			}
			System.out.println("dotp_cluster=" + dotp_clus);
			mean_projection[0] = 0;
			mean_projection[1] = 0;
			mean_projection[2] = 0;
			mean_projection[3] = 0;
			mean_projection[4] = 0;

			for (int i = 0; i < vocabulary.size(); i++) {
				mean_projection[0] += cluster1[0][i] * mean[i];
				mean_projection[1] += cluster1[1][i] * mean[i];
				mean_projection[2] += cluster1[2][i] * mean[i];
				mean_projection[3] += cluster1[3][i] * mean[i];
				mean_projection[4] += cluster1[4][i] * mean[i];
			}
			System.out.println(mean_projection[0]);
			System.out.println(mean_projection[1]);
			System.out.println(mean_projection[2]);
			int[] cluster_size = new int[5];
			TObjectIntHashMap<String>[] class_clusters = new TObjectIntHashMap[5];
			for (int i = 0; i < 5; i++) {
				class_clusters[i] = new TObjectIntHashMap<String>();
			}
			TreeSet<DescSorter> sortedRel1 = new TreeSet<DescSorter>();
			TreeSet<DescSorter> sortedRel2 = new TreeSet<DescSorter>();
			// ArrayList<TreeSet<DescSorter>> sortedRel = new
			// ArrayList<TreeSet<DescSorter>>(
			// 2);
			// for (int i = 0; i < 2; i++)
			// sortedRel.add(i, new TreeSet<DescSorter>());
			TObjectIntHashMap<String> wrdC = new TObjectIntHashMap<String>();
			double[][] dotps = new double[3][this.docs.size()];

			// Dotps
			for (int j = 0; j < this.docs.size(); j++) {

				double dotp = 0;

				for (int i = 0; i < 3; i++) {
					int[] keys = this.docs.get(j).normed_tfVector.keys();
					dotp = 0;
					for (int k = 0; k < keys.length; k++) {
						// if(i==0)
						dotp += this.docs.get(j).normed_tfVector.get(keys[k])
								* this.normalizedEignVector[i][keys[k]];
						// if(i==1)
						// dotp += this.docs.get(j).normed_tfVector.get(keys[k])
						// * cluster1[0][keys[k]];
						// dotp += this.docs.get(j).normed_tfVector.get(keys[k])
						// * cluster1[1][keys[k]];
					}
					dotp -= this.mean_projection[i];
					dotps[i][j] = dotp;
				}
			}
			BufferedWriter dot_products = null;
			try {
				System.out.println(" Writing");
				dot_products = new BufferedWriter(new FileWriter("./dotps"));
				for (int j = 0; j < this.docs.size(); j++) {
					String name = this.docs.get(j).docName;
					dot_products.write(name + ":");
					for (int i = 0; i < 3; i++) {
						dot_products.write(dotps[i][j] + ":");
					}
					dot_products.write("\n");
				}

				dot_products.close();
				System.out.println(" Writing Done");

			} catch (Exception e) {
				e.printStackTrace();
				System.exit(0);
			}
			// end

			for (int j = 0; j < this.docs.size(); j++) {
				int max_cluster = -1;
				double dotp = 0;
				double max_dotp = Double.NEGATIVE_INFINITY;
				for (int i = 1; i < 4; i++) {
					int[] keys = this.docs.get(j).normed_tfVector.keys();
					dotp = 0;
					for (int k = 0; k < keys.length; k++) {
						// if(i==0)
						dotp += this.docs.get(j).normed_tfVector.get(keys[k])
								* cluster1[i][keys[k]];
						// if(i==1)
						// dotp += this.docs.get(j).normed_tfVector.get(keys[k])
						// * cluster1[0][keys[k]];
						// dotp += this.docs.get(j).normed_tfVector.get(keys[k])
						// * cluster1[1][keys[k]];
					}
					dotp -= this.mean_projection[i];
					// dotp -= this.mean_projection[1];
					// dotp *= dotp;
					String name = this.docs.get(j).docName;
					// System.out.println(name+":"+dotp);
					String labl = name.split("/")[0];
					// sortedRel.get(i).add(new DescSorter((String) labl,
					// dotp));
					/*
					 * if (dotp < 0) { max_cluster = 0; } else max_cluster = 1;
					 */

					if (dotp >= max_dotp) {
						max_cluster = i;
						max_dotp = dotp;
					}

				}

				if (max_cluster != -1) {
					String name = this.docs.get(j).docName;
					wrdC.put(name.toString(),
							this.docs.get(j).word_counts.size());
					if (max_cluster == 0) {
						sortedRel1.add(new DescSorter((String) name, dotp));
					}
					if (max_cluster == 1) {
						sortedRel2.add(new DescSorter((String) name, dotp));
					}
					String labl = name.split("/")[0];
					class_clusters[max_cluster].adjustOrPutValue(labl, 1, 1);
					cluster_size[max_cluster]++;

				} else {
					// System.out.println("Error in finding relevance cluster!!");
				}

			}
			Iterator<DescSorter> DescIter = sortedRel1.iterator();
			while (DescIter.hasNext()) {
				DescSorter info = DescIter.next();

				System.out.println("Doc=" + info.getID() + " :: dotp="
						+ info.getWeight() + " :: length="
						+ wrdC.get(info.getID()));
			}
			System.out.println("-------------------------------------------");
			DescIter = sortedRel2.iterator();
			while (DescIter.hasNext()) {
				DescSorter info = DescIter.next();
				System.out.println("Doc=" + info.getID() + " :: dotp="
						+ info.getWeight() + " :: length="
						+ wrdC.get(info.getID()));
			}
			/*
			 * TObjectIntHashMap<String>[] topRel = new TObjectIntHashMap[2];
			 * for (int i = 0; i < 2; i++) topRel[i] = new
			 * TObjectIntHashMap<String>(); Iterator<DescSorter> DescIter =
			 * sortedRel.get(0).iterator(); for (int i = 0; i < 100; i++) {
			 * DescSorter info = DescIter.next();
			 * topRel[0].adjustOrPutValue(info.getID(), 1, 1); } DescIter =
			 * sortedRel.get(1).iterator(); for (int i = 0; i < 100; i++) {
			 * DescSorter info = DescIter.next();
			 * topRel[1].adjustOrPutValue(info.getID(), 1, 1); }
			 */

			output.write("\n\n-----------Document Distribution in clusters--------------\n");
			for (int i = 0; i < num_clusters; i++) {

				TreeSet<DescSorter> sortedCat = new TreeSet<DescSorter>();
				Object[] keys = class_clusters[i].keys();
				for (Object key : keys) {
					sortedCat.add(new DescSorter((String) key,
							class_clusters[i].get((String) key)));
				}
				DescIter = sortedCat.iterator();
				output.write("Cluster-" + i + ": # docs in cluster= "
						+ cluster_size[i] + ":\n");
				while (DescIter.hasNext()) {
					DescSorter info = DescIter.next();
					output.write(info.getID() + ":" + info.getWeight() + "\t");
				}
				output.write("\n\n");
			}
			/*
			 * output.write(
			 * "\n\n-----------Top documents in clusters--------------\n"); for
			 * (int i = 0; i < num_clusters; i++) {
			 * 
			 * output.write("Cluster-" + i + ":\n"); Object[] keys =
			 * topRel[i].keys(); for (int j = 0; j < keys.length; j++) {
			 * output.write((String) keys[j] + ":" + topRel[i].get((String)
			 * keys[j]) + "\t"); }
			 * 
			 * output.write("\n\n"); }
			 */
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

	}

	public void printClusterDocRelevance(BufferedWriter output,
			int num_clusters, int count) {
		try {
			for (int i = 0; i < num_clusters; i++) {
				TreeSet<IDSorter> sortedWords = new TreeSet<IDSorter>();
				for (int j = 0; j < this.docs.size(); j++) {
					int[] keys = this.docs.get(j).normed_tfidfVector.keys();
					double dotp = 0;
					for (int k = 0; k < keys.length; k++) {
						dotp += this.docs.get(j).normed_tfidfVector
								.get(keys[k])
								* this.cluster_states[i].get(keys[k]);
					}

					sortedWords.add(new IDSorter(j, dotp));

				}
				TObjectIntHashMap<String> lablCount = new TObjectIntHashMap<String>();
				Iterator<IDSorter> iterator = sortedWords.iterator();

				int docs = 1;
				// while(iterator.hasNext() && words<keys.length-count){
				// IDSorter info = iterator.next();
				// words++;
				// }
				output.write("\n------------\n");
				output.write("Cluster-" + i + ": Prob:" + this.cluster_probs[i]
						+ "\n");
				while (iterator.hasNext() && docs < 500) {
					IDSorter info = iterator.next();
					String name = this.docs.get(info.getID()).docName;
					String labl = name.split("/")[0];
					lablCount.adjustOrPutValue(labl, 1, 1);
					// if (docs < count)
					// output.write(this.docs.get(info.getID()).docName + ":"
					// + info.getWeight() + "\t");
					docs++;

				}
				output.write("\n");
				Object[] keys = lablCount.keys();
				TreeSet<DescSorter> sortedCat = new TreeSet<DescSorter>();

				for (Object key : keys) {
					sortedCat.add(new DescSorter((String) key, lablCount
							.get((String) key)));
				}
				Iterator<DescSorter> DescIter = sortedCat.iterator();
				while (DescIter.hasNext()) {
					DescSorter info = DescIter.next();
					output.write(info.getID() + ":" + info.getWeight() + "\t");
				}
				// System.out.println("-----"+words+"-------");

			}
			int[] cluster_size = new int[num_clusters];
			TObjectIntHashMap<String>[] class_clusters = new TObjectIntHashMap[num_clusters];
			for (int i = 0; i < num_clusters; i++) {
				class_clusters[i] = new TObjectIntHashMap<String>();
			}

			for (int j = 0; j < this.docs.size(); j++) {
				int max_cluster = -1;
				double max_dotp = 0;
				for (int i = 0; i < num_clusters; i++) {
					int[] keys = this.docs.get(j).normed_tfidfVector.keys();
					double dotp = 0;
					for (int k = 0; k < keys.length; k++) {
						dotp += this.docs.get(j).normed_tfidfVector
								.get(keys[k])
								* this.cluster_states[i].get(keys[k]);
					}
					if (dotp > max_dotp) {
						max_cluster = i;
						max_dotp = dotp;
					}

				}
				if (max_cluster != -1) {
					String name = this.docs.get(j).docName;
					String labl = name.split("/")[0];
					class_clusters[max_cluster].adjustOrPutValue(labl, 1, 1);
					cluster_size[max_cluster]++;
				} else {
					// System.out.println("Error in finding relevance cluster!!");
				}

			}
			output.write("\n\n");
			for (int i = 0; i < num_clusters; i++) {

				TreeSet<DescSorter> sortedCat = new TreeSet<DescSorter>();
				Object[] keys = class_clusters[i].keys();
				output.write("Cluster-" + i + ": Prob:" + this.cluster_probs[i]
						+ ": # docs in cluster= " + cluster_size[i]
						+ "\n class distribution: ");
				for (Object key : keys) {
					output.write(key.toString() + ":"
							+ class_clusters[i].get((String) key) + "\t");
				}

				output.write("\n\n");
			}

			output.write("\n------------\n");

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

	}

	public void calculateEntropy(BufferedWriter output, int num_clusters) {
		try {
			int[] cluster_size = new int[num_clusters];
			TObjectIntHashMap<String>[] class_clusters = new TObjectIntHashMap[num_clusters];
			TObjectIntHashMap<String>[] super_class_clusters = new TObjectIntHashMap[num_clusters];
			for (int i = 0; i < num_clusters; i++) {
				class_clusters[i] = new TObjectIntHashMap<String>();
				super_class_clusters[i] = new TObjectIntHashMap<String>();
			}

			for (int j = 0; j < this.docs.size(); j++) {
				int max_cluster = -1;
				double max_dotp = 0;
				for (int i = 0; i < num_clusters; i++) {
					int[] keys = this.docs.get(j).normed_tfidfVector.keys();
					double dotp = 0;
					for (int k = 0; k < keys.length; k++) {
						dotp += this.docs.get(j).normed_tfidfVector
								.get(keys[k])
								* this.cluster_states[i].get(keys[k]);
					}
					if (dotp > max_dotp) {
						max_cluster = i;
						max_dotp = dotp;
					}

				}
				if (max_cluster != -1) {
					String name = this.docs.get(j).docName;
					String labl = name.split("/")[0];
					class_clusters[max_cluster].adjustOrPutValue(labl, 1, 1);
					String super_labl = labl.split("\\.")[0];
					super_class_clusters[max_cluster].adjustOrPutValue(
							super_labl, 1, 1);
					cluster_size[max_cluster]++;
				} else {
					// System.out.println("Error in finding relevance cluster!!");
				}

			}
			output.write("\n\n-----------Document Distribution in clusters--------------\n");
			for (int i = 0; i < num_clusters; i++) {

				TreeSet<DescSorter> sortedCat = new TreeSet<DescSorter>();
				Object[] keys = class_clusters[i].keys();
				for (Object key : keys) {
					sortedCat.add(new DescSorter((String) key,
							class_clusters[i].get((String) key)));
				}
				Iterator<DescSorter> DescIter = sortedCat.iterator();
				output.write("Cluster-" + i + ": Prob:" + this.cluster_probs[i]
						+ ": # docs in cluster= " + cluster_size[i] + ":\n");
				while (DescIter.hasNext()) {
					DescSorter info = DescIter.next();
					output.write(info.getID() + ":" + info.getWeight() + "\t");
				}
				output.write("\n\n");
			}
			output.write("\n\n-----------Cluster Entropy for sub class--------------\n");
			double entropy = 0;
			for (int i = 0; i < num_clusters; i++) {
				double cluster_entropy = 0;
				if (cluster_size[i] > 0) {
					for (Object key : class_clusters[i].keys()) {
						cluster_entropy += (class_clusters[i].get((String) key) * Math
								.log((double) class_clusters[i]
										.get((String) key)
										/ (double) cluster_size[i]));
						// System.out.println(class_clusters[i].get((String)
						// key)+":"+cluster_size[i]);
					}

					cluster_entropy = (-1) * cluster_entropy
							/ (double) this.docs.size();

				}
				output.write("Cluster-" + i + ": Prob:" + this.cluster_probs[i]
						+ ": # docs in cluster= " + cluster_size[i] + ":");

				output.write("Class Entropy :" + cluster_entropy);

				output.write("\n\n");
				entropy += cluster_entropy;
			}
			output.write("Overall Entropy :" + entropy);

			output.write("\n\n-----------Cluster Entropy for super class--------------\n");
			entropy = 0;
			for (int i = 0; i < num_clusters; i++) {
				double cluster_entropy = 0;
				if (cluster_size[i] > 0) {
					for (Object key : super_class_clusters[i].keys()) {
						cluster_entropy += (super_class_clusters[i]
								.get((String) key) * Math
								.log((double) super_class_clusters[i]
										.get((String) key)
										/ (double) cluster_size[i]));
					}
					cluster_entropy = (-1) * cluster_entropy
							/ (double) this.docs.size();

				}
				output.write("Cluster-" + i + ": Prob:" + this.cluster_probs[i]
						+ ": # docs in cluster= " + cluster_size[i] + ":");

				output.write("Class Entropy :" + cluster_entropy);

				output.write("\n\n");
				entropy += cluster_entropy;
			}
			output.write("Overall Entropy :" + entropy + "\n");
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

	}

	public void calculateEntropyHardSubRelevance(BufferedWriter output,
			int num_clusters) {
		try {
			int[] cluster_size = new int[num_clusters];
			TObjectIntHashMap<String>[] class_clusters = new TObjectIntHashMap[num_clusters];
			TObjectIntHashMap<String>[] super_class_clusters = new TObjectIntHashMap[num_clusters];
			for (int i = 0; i < num_clusters; i++) {
				class_clusters[i] = new TObjectIntHashMap<String>();
				super_class_clusters[i] = new TObjectIntHashMap<String>();
			}

			for (int j = 0; j < this.docs.size(); j++) {
				int max_cluster = -1;
				double max_dotp = 0;
				for (int i = 0; i < 2; i++) {
					int[] keys = this.docs.get(j).normed_tfidfVector.keys();
					double dotp = 0;
					for (int k = 0; k < keys.length; k++) {
						dotp += this.docs.get(j).normed_tfidfVector
								.get(keys[k])
								* this.subspace_vectors[2 * this.doc_cluster[j]
										+ i].get(keys[k]);
					}
					if (dotp * dotp > max_dotp) {
						max_cluster = i;
						max_dotp = dotp * dotp;
					}

				}
				if (max_cluster != -1) {
					String name = this.docs.get(j).docName;
					String labl = name.split("/")[0];
					class_clusters[2 * this.doc_cluster[j] + max_cluster]
							.adjustOrPutValue(labl, 1, 1);
					String super_labl = labl.split("\\.")[0];
					super_class_clusters[2 * this.doc_cluster[j] + max_cluster]
							.adjustOrPutValue(super_labl, 1, 1);
					cluster_size[2 * this.doc_cluster[j] + max_cluster]++;
				} else {
					// System.out.println("Error in finding relevance cluster!!");
				}

			}
			output.write("\n\n-----------Document Distribution in clusters--------------\n");
			for (int i = 0; i < num_clusters; i++) {

				TreeSet<DescSorter> sortedCat = new TreeSet<DescSorter>();
				Object[] keys = class_clusters[i].keys();
				for (Object key : keys) {
					sortedCat.add(new DescSorter((String) key,
							class_clusters[i].get((String) key)));
				}
				Iterator<DescSorter> DescIter = sortedCat.iterator();
				output.write("Cluster-" + i + ": Prob:" + this.alphas[i]
						+ ": # docs in cluster= " + cluster_size[i] + ":\n");
				while (DescIter.hasNext()) {
					DescSorter info = DescIter.next();
					output.write(info.getID() + ":" + info.getWeight() + "\t");
				}
				output.write("\n\n");
			}
			/*
			 * output.write(
			 * "\n\n-----------Cluster Entropy for sub class--------------\n");
			 * double entropy = 0; for (int i = 0; i < num_clusters; i++) {
			 * double cluster_entropy = 0; if (cluster_size[i] > 0) { for
			 * (Object key : class_clusters[i].keys()) { cluster_entropy +=
			 * (class_clusters[i].get((String) key) * Math .log((double)
			 * class_clusters[i] .get((String) key) / (double)
			 * cluster_size[i])); //
			 * System.out.println(class_clusters[i].get((String) //
			 * key)+":"+cluster_size[i]); }
			 * 
			 * cluster_entropy = (-1) * cluster_entropy / (double)
			 * this.docs.size();
			 * 
			 * } output.write("Cluster-" + i + ": Prob:" + this.cluster_probs[i]
			 * + ": # docs in cluster= " + cluster_size[i] + ":");
			 * 
			 * output.write("Class Entropy :" + cluster_entropy);
			 * 
			 * output.write("\n\n"); entropy += cluster_entropy; }
			 * output.write("Overall Entropy :" + entropy);
			 * 
			 * output.write(
			 * "\n\n-----------Cluster Entropy for super class--------------\n"
			 * ); entropy = 0; for (int i = 0; i < num_clusters; i++) { double
			 * cluster_entropy = 0; if (cluster_size[i] > 0) { for (Object key :
			 * super_class_clusters[i].keys()) { cluster_entropy +=
			 * (super_class_clusters[i] .get((String) key) * Math .log((double)
			 * super_class_clusters[i] .get((String) key) / (double)
			 * cluster_size[i])); } cluster_entropy = (-1) * cluster_entropy /
			 * (double) this.docs.size();
			 * 
			 * } output.write("Cluster-" + i + ": Prob:" + this.cluster_probs[i]
			 * + ": # docs in cluster= " + cluster_size[i] + ":");
			 * 
			 * output.write("Class Entropy :" + cluster_entropy);
			 * 
			 * output.write("\n\n"); entropy += cluster_entropy; }
			 * output.write("Overall Entropy :" + entropy + "\n");
			 */
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

	}

	public void calculateEntropyHardRelevance(BufferedWriter output,
			int num_clusters) {
		try {
			int[] cluster_size = new int[num_clusters];
			TObjectIntHashMap<String>[] class_clusters = new TObjectIntHashMap[num_clusters];
			TObjectIntHashMap<String>[] super_class_clusters = new TObjectIntHashMap[num_clusters];
			for (int i = 0; i < num_clusters; i++) {
				class_clusters[i] = new TObjectIntHashMap<String>();
				super_class_clusters[i] = new TObjectIntHashMap<String>();
			}

			for (int j = 0; j < this.docs.size(); j++) {
				int max_cluster = -1;
				double max_dotp = 0;
				for (int i = 0; i < num_clusters; i++) {
					int[] keys = this.docs.get(j).normed_tfidfVector.keys();
					double dotp = 0;
					for (int k = 0; k < keys.length; k++) {
						dotp += this.docs.get(j).normed_tfidfVector
								.get(keys[k])
								* this.cluster_states[i].get(keys[k]);
					}
					if (dotp > max_dotp) {
						max_cluster = i;
						max_dotp = dotp;
					}

				}
				if (max_cluster != -1) {
					String name = this.docs.get(j).docName;
					String labl = name.split("/")[0];
					class_clusters[max_cluster].adjustOrPutValue(labl, 1, 1);
					String super_labl = labl.split("\\.")[0];
					super_class_clusters[max_cluster].adjustOrPutValue(
							super_labl, 1, 1);
					cluster_size[max_cluster]++;
				} else {
					// System.out.println("Error in finding relevance cluster!!");
				}

			}
			output.write("\n\n-----------Document Distribution in clusters--------------\n");
			for (int i = 0; i < num_clusters; i++) {

				TreeSet<DescSorter> sortedCat = new TreeSet<DescSorter>();
				Object[] keys = class_clusters[i].keys();
				for (Object key : keys) {
					sortedCat.add(new DescSorter((String) key,
							class_clusters[i].get((String) key)));
				}
				Iterator<DescSorter> DescIter = sortedCat.iterator();
				output.write("Cluster-" + i + ": Prob:" + this.cluster_probs[i]
						+ ": # docs in cluster= " + cluster_size[i] + ":\n");
				while (DescIter.hasNext()) {
					DescSorter info = DescIter.next();
					output.write(info.getID() + ":" + info.getWeight() + "\t");
				}
				output.write("\n\n");
			}
			output.write("\n\n-----------Cluster Entropy for sub class--------------\n");
			double entropy = 0;
			for (int i = 0; i < num_clusters; i++) {
				double cluster_entropy = 0;
				if (cluster_size[i] > 0) {
					for (Object key : class_clusters[i].keys()) {
						cluster_entropy += (class_clusters[i].get((String) key) * Math
								.log((double) class_clusters[i]
										.get((String) key)
										/ (double) cluster_size[i]));
						// System.out.println(class_clusters[i].get((String)
						// key)+":"+cluster_size[i]);
					}

					cluster_entropy = (-1) * cluster_entropy
							/ (double) this.docs.size();

				}
				output.write("Cluster-" + i + ": Prob:" + this.cluster_probs[i]
						+ ": # docs in cluster= " + cluster_size[i] + ":");

				output.write("Class Entropy :" + cluster_entropy);

				output.write("\n\n");
				entropy += cluster_entropy;
			}
			output.write("Overall Entropy :" + entropy);

			output.write("\n\n-----------Cluster Entropy for super class--------------\n");
			entropy = 0;
			for (int i = 0; i < num_clusters; i++) {
				double cluster_entropy = 0;
				if (cluster_size[i] > 0) {
					for (Object key : super_class_clusters[i].keys()) {
						cluster_entropy += (super_class_clusters[i]
								.get((String) key) * Math
								.log((double) super_class_clusters[i]
										.get((String) key)
										/ (double) cluster_size[i]));
					}
					cluster_entropy = (-1) * cluster_entropy
							/ (double) this.docs.size();

				}
				output.write("Cluster-" + i + ": Prob:" + this.cluster_probs[i]
						+ ": # docs in cluster= " + cluster_size[i] + ":");

				output.write("Class Entropy :" + cluster_entropy);

				output.write("\n\n");
				entropy += cluster_entropy;
			}
			output.write("Overall Entropy :" + entropy + "\n");
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

	}

	public void calculateEntropyKMeans(BufferedWriter output, int num_clusters) {
		try {
			int[] cluster_size = new int[num_clusters];
			TObjectIntHashMap<String>[] class_clusters = new TObjectIntHashMap[num_clusters];
			TObjectIntHashMap<String>[] super_class_clusters = new TObjectIntHashMap[num_clusters];
			for (int i = 0; i < num_clusters; i++) {
				class_clusters[i] = new TObjectIntHashMap<String>();
				super_class_clusters[i] = new TObjectIntHashMap<String>();
			}
			double[] cluster_norm = new double[num_clusters];
			for (int j = 0; j < num_clusters; j++) {
				for (int v = 0; v < num_clusters; v++) {
					cluster_norm[j] += Math.pow(this.cluster_states[j].get(v),
							2);
				}
			}

			for (int i = 0; i < this.docs.size(); i++) {
				double min = Double.MAX_VALUE;
				int max_cluster = -1;
				int[] keys_data = this.docs.get(i).normed_tfidfVector.keys();
				for (int l = 0; l < num_clusters; l++) {
					double dist = cluster_norm[l];
					for (int k = 0; k < keys_data.length; k++) {
						dist -= Math.pow(
								this.cluster_states[l].get(keys_data[k]), 2);
						dist += Math
								.pow((this.cluster_states[l].get(keys_data[k]) - this.docs
										.get(i).normed_tfidfVector
										.get(keys_data[k])), 2);
					}
					if (dist < min) {
						min = dist;
						max_cluster = l;
					}
				}

				String name = this.docs.get(i).docName;
				String labl = name.split("/")[0];
				class_clusters[max_cluster].adjustOrPutValue(labl, 1, 1);
				String super_labl = labl.split("\\.")[0];
				super_class_clusters[max_cluster].adjustOrPutValue(super_labl,
						1, 1);
				cluster_size[max_cluster]++;

				// doc_cluster[i] = max_cluster;
				// cluster_size[max_cluster]++;

			}

			output.write("\n\n-----------Document Distribution in clusters--------------\n");
			for (int i = 0; i < num_clusters; i++) {

				TreeSet<DescSorter> sortedCat = new TreeSet<DescSorter>();
				Object[] keys = class_clusters[i].keys();
				for (Object key : keys) {
					sortedCat.add(new DescSorter((String) key,
							class_clusters[i].get((String) key)));
				}
				Iterator<DescSorter> DescIter = sortedCat.iterator();
				output.write("Cluster-" + i + ": # docs in cluster= "
						+ cluster_size[i] + ":\n");
				while (DescIter.hasNext()) {
					DescSorter info = DescIter.next();
					output.write(info.getID() + ":" + info.getWeight() + "\t");
				}
				output.write("\n\n");
			}
			output.write("\n\n-----------Cluster Entropy for sub class--------------\n");
			double entropy = 0;
			for (int i = 0; i < num_clusters; i++) {
				double cluster_entropy = 0;
				if (cluster_size[i] > 0) {
					for (Object key : class_clusters[i].keys()) {
						cluster_entropy += (class_clusters[i].get((String) key) * Math
								.log((double) class_clusters[i]
										.get((String) key)
										/ (double) cluster_size[i]));
						// System.out.println(class_clusters[i].get((String)
						// key)+":"+cluster_size[i]);
					}

					cluster_entropy = (-1) * cluster_entropy
							/ (double) this.docs.size();

				}
				output.write("Cluster-" + i + ": # docs in cluster= "
						+ cluster_size[i] + ":");

				output.write("Class Entropy :" + cluster_entropy);

				output.write("\n\n");
				entropy += cluster_entropy;
			}
			output.write("Overall Entropy :" + entropy);

			output.write("\n\n-----------Cluster Entropy for super class--------------\n");
			entropy = 0;
			for (int i = 0; i < num_clusters; i++) {
				double cluster_entropy = 0;
				if (cluster_size[i] > 0) {
					for (Object key : super_class_clusters[i].keys()) {
						cluster_entropy += (super_class_clusters[i]
								.get((String) key) * Math
								.log((double) super_class_clusters[i]
										.get((String) key)
										/ (double) cluster_size[i]));
					}
					cluster_entropy = (-1) * cluster_entropy
							/ (double) this.docs.size();

				}
				output.write("Cluster-" + i + ": # docs in cluster= "
						+ cluster_size[i] + ":");

				output.write("Class Entropy :" + cluster_entropy);

				output.write("\n\n");
				entropy += cluster_entropy;
			}
			output.write("Overall Entropy :" + entropy + "\n");
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

	}

	public void calculateEntropyRandom(BufferedWriter output, int num_clusters) {
		try {

			double[] cluster_ent = new double[num_clusters];
			double[] super_cluster_ent = new double[num_clusters];
			int count = 0;
			for (int iters = 1; iters < 50; iters++) {
				count++;
				Random r = new Random(iters);
				int[] cluster_size = new int[num_clusters];
				TObjectIntHashMap<String>[] class_clusters = new TObjectIntHashMap[num_clusters];
				TObjectIntHashMap<String>[] super_class_clusters = new TObjectIntHashMap[num_clusters];
				for (int i = 0; i < num_clusters; i++) {
					class_clusters[i] = new TObjectIntHashMap<String>();
					super_class_clusters[i] = new TObjectIntHashMap<String>();
				}
				for (int i = 0; i < this.docs.size(); i++) {
					int max_cluster = r.nextInt(num_clusters);
					String name = this.docs.get(i).docName;
					String labl = name.split("/")[0];
					class_clusters[max_cluster].adjustOrPutValue(labl, 1, 1);
					String super_labl = labl.split("\\.")[0];
					super_class_clusters[max_cluster].adjustOrPutValue(
							super_labl, 1, 1);
					cluster_size[max_cluster]++;
				}
				for (int i = 0; i < num_clusters; i++) {
					double cluster_entropy = 0;
					if (cluster_size[i] > 0) {
						for (Object key : class_clusters[i].keys()) {
							cluster_entropy += (class_clusters[i]
									.get((String) key) * Math
									.log((double) class_clusters[i]
											.get((String) key)
											/ (double) cluster_size[i]));
						}
						cluster_entropy = (-1) * cluster_entropy
								/ (double) this.docs.size();
					}
					cluster_ent[i] += cluster_entropy;
				}
				for (int i = 0; i < num_clusters; i++) {
					double cluster_entropy = 0;
					if (cluster_size[i] > 0) {
						for (Object key : super_class_clusters[i].keys()) {
							cluster_entropy += (super_class_clusters[i]
									.get((String) key) * Math
									.log((double) super_class_clusters[i]
											.get((String) key)
											/ (double) cluster_size[i]));
						}
						cluster_entropy = (-1) * cluster_entropy
								/ (double) this.docs.size();
					}
					super_cluster_ent[i] += cluster_entropy;
				}
			}
			double ent = 0, sup_ent = 0;
			for (int i = 0; i < num_clusters; i++) {
				cluster_ent[i] /= count;
				super_cluster_ent[i] /= count;
				ent += cluster_ent[i];
				sup_ent += super_cluster_ent[i];
			}
			output.write("\n\n-----------Cluster Entropy for sub class--------------\n");
			for (int i = 0; i < num_clusters; i++) {
				output.write("Cluster-" + i + ":");
				output.write("Class Entropy :" + cluster_ent[i]);
				output.write("\n\n");
			}
			output.write("Overall Entropy :" + ent);
			output.write("\n\n-----------Cluster Entropy for super class--------------\n");
			for (int i = 0; i < num_clusters; i++) {
				output.write("Cluster-" + i + ":");
				output.write("Class Entropy :" + super_cluster_ent[i]);
				output.write("\n\n");
			}
			output.write("Overall Entropy :" + sup_ent + "\n");
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}

	}

	public double reconstruction_error(int num_clusters, double[] doc_proj_sum,
			double[] factor_p_sum) {
		double ll = 0.0, dotp, cum_resp;
		// for(int j=0;j<num_clusters;j++)
		// System.out.println(this.cluster_probs[j]);
		for (int i = 0; i < this.docs.size(); i++) {

			cum_resp = 0;
			int[] keys_data = this.docs.get(i).normed_tfidfVector.keys();
			for (int j = 0; j < num_clusters; j++) {
				dotp = 0;
				for (int k = 0; k < keys_data.length; k++) {
					if (this.cluster_states[j].containsKey(keys_data[k]))
						dotp += this.docs.get(i).normed_tfidfVector
								.get(keys_data[k])
								* this.cluster_states[j].get(keys_data[k]);
				}
				cum_resp += dotp * dotp * this.cluster_probs[j];

			}
			if (cum_resp > 0)
				ll += Math.log(cum_resp);
		}
		return ll;

	}

	public void saveState(String file) {
		try {
			FileOutputStream fout = new FileOutputStream(file);
			ObjectOutputStream oos = new ObjectOutputStream(fout);
			this.writeObject(oos);
			oos.close();
			fout.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public void read(String f) {

		try {
			FileInputStream fin = new FileInputStream(f);
			ObjectInputStream ois = new ObjectInputStream(fin);
			this.readObject(ois);
			// this.initializeTypeTopicCounts(); // To work around a bug in
			// Trove?
			ois.close();
		} catch (IOException e) {
			System.err.println("Exception reading file " + f + ": " + e);
		} catch (ClassNotFoundException e) {
			System.err.println("Exception reading file " + f + ": " + e);
		}

	}

	private static final long serialVersionUID = 1;
	private static final int CURRENT_SERIAL_VERSION = 0;

	private void writeObject(ObjectOutputStream out) throws IOException {
		out.writeInt(CURRENT_SERIAL_VERSION);
		out.writeObject(vocabulary);
		// out.writeInt(size);
		vocab_size = vocabulary.size();
		out.writeInt(vocab_size);
		out.writeObject(idf);
		for (int i = 0; i < docs.size(); i++)
			out.writeObject(docs.get(i));
	}

	private void readObject(ObjectInputStream in) throws IOException,
			ClassNotFoundException {
		int version = in.readInt();
		vocabulary = (Alphabet) in.readObject();
		// size = in.readInt();
		vocab_size = in.readInt();
		idf = (TIntDoubleHashMap) in.readObject();
		docs = new Vector<Document>();

		// for (int i = 0; i < size; i++)
		// docs.add((Document) in.readObject());

	}

	// matrix multiplication routine for arpack

	/*
	 * int getLocalTopicCounts(int word, int topic){
	 * 
	 * } int getGlobalTopicCounts(int topic){
	 * 
	 * }
	 */

}
