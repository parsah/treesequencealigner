package model;

import java.util.Set;
import java.util.TreeMap;

/**
 * A NodeMap references mappings between various nodes in cases where
 * a key aligns perfectly to a value.
 * */
public class NodeType extends TreeMap<String, String> {
	private static final long serialVersionUID = 1L;

	public NodeType() {
		this.put("A", "A");
		this.put("T", "T");
		this.put("C", "C");
	}
	
	/**
	 * Returns the keys as-part of this nodemap collection.
	 * @return Set key-collection comprising alignment nodes.
	 * */
	public Set<String> keys() {
		return this.keySet();
	}
	
	/**
	 * Return number of entries in the NodeMap
	 * @return int references number of keys in the map.
	 * */
	public int getSize() {
		return this.getSize();
	}

}
