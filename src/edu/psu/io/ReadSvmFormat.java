package edu.psu.io;

import java.io.*;

public class ReadSvmFormat implements Serializable {
	
	// Serialization

	private static final long serialVersionUID = 1;
	private static final int CURRENT_SERIAL_VERSION = 1;

	private void writeObject (ObjectOutputStream out) throws IOException {
		/*out.writeInt (CURRENT_SERIAL_VERSION);
		out.writeInt (entries.size());
		for (int i = 0; i < entries.size(); i++)
			out.writeObject (entries.get(i));
		out.writeBoolean (growthStopped);
		out.writeObject (entryClass);
		out.writeObject(instanceId);*/
	}

	private void readObject (ObjectInputStream in) throws IOException, ClassNotFoundException {
		/*int version = in.readInt ();
		int size = in.readInt();
		entries = new ArrayList (size);
		map = new gnu.trove.TObjectIntHashMap (size);
		for (int i = 0; i < size; i++) {
			Object o = in.readObject();
			map.put (o, i);
			entries. add (o);
		}
		growthStopped = in.readBoolean();
		entryClass = (Class) in.readObject();
		if (version >0 ){ // instanced id added in version 1S
			instanceId = (VMID) in.readObject();
		}*/
	}

	public static void main(String[] args) {
		//To-Do
		
		
		
	}

}
