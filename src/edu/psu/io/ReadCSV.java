package edu.psu.io;

import java.io.*;
import edu.psu.types.*;
import java.util.*;

public class ReadCSV {
	
	// Serialization

	private static final long serialVersionUID = 1;
	private static final int CURRENT_SERIAL_VERSION = 1;
	//read from a file containing document per line ---- needs changes in the document constructor
	public static Vector<Document> readSequence(String file, boolean isLable) {
		Vector<Document> docs = new Vector<Document>();
		String line = null;
		try {
			BufferedReader br = null;
			br = new BufferedReader(new FileReader(file));
			while ((line = br.readLine()) != null) {
				Document doc = new Document();
				docs.add(doc);
			}
			
		} catch (java.io.IOException e) {
			throw new IllegalArgumentException("IOException " + e);
		}

		return docs;

	}
	
	//read from a file containing document per line ---- needs changes in the document constructor
	public static Vector<Document> readDir(String file, boolean isLable) {
		Vector<Document> docs = new Vector<Document>();
		String line = null;
		try {
			BufferedReader br = null;
			br = new BufferedReader(new FileReader(file));
			while ((line = br.readLine()) != null) {
				Document doc = new Document();
				docs.add(doc);
			}
			
		} catch (java.io.IOException e) {
			throw new IllegalArgumentException("IOException " + e);
		}

		return docs;

	}

	public static Vector<Document> readSequence(String file) {

		return readSequence(file, true);

	}

	static Vector<Document> readBag(String file) {
		Vector<Document> docs = null;
		return docs;
	}

	public static void main(String[] args) {
		
		if (args.length == 0) {
			System.out.println("USAGE: ReadCSV inputfile output_model_file");
		}
		String input = args[0].toString();
		String output = args[1].toString();
		Corpus corpus = new Corpus();
		corpus.addDocuments(readSequence(input));
		corpus.normalizeCorpus();
		corpus.saveState(output);
		Corpus corpus1 = new Corpus(); 
		corpus1.read(output);
		System.out.println(corpus1.vocab_size);
		System.out.println(corpus1.docs.size());
		System.out.println(corpus1.docs.size());
		Document doc = corpus1.docs.get(0);
		System.out.println(doc.docName);
		
		
	}

}
