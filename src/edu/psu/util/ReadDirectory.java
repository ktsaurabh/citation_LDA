package edu.psu.util;

import java.util.*;
import java.io.*;

public class ReadDirectory {
	public static String[] list(String path) {
		File dir = new File(path);
		String absPath = dir.getAbsolutePath();
		FilenameFilter filter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return !name.startsWith(".");
			}
		};
		String[] children = dir.list(filter);

		if (children == null) {
			return null;
		} else {
			for (int i = 0; i < children.length; i++) {
				// Get filename of file or directory
				children[i]=absPath+"/"+children[i];
			}
		}
		return children;

		// It is also possible to filter the list of returned files.
		// This example does not return any files that start with `.'.

	}
	public static void main(String[] args){
		ReadDirectory.list(args[0]);
		
	}
	
}
