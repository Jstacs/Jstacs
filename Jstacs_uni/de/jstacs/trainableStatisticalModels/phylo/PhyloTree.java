package de.jstacs.trainableStatisticalModels.phylo;

import de.jstacs.NonParsableException;
import de.jstacs.Storable;
import de.jstacs.io.XMLParser;
import java.util.ArrayList;

/**
 * This class implements a simple (phylogenetic) tree.
 * A PhyloTree contains the root node of type {@link PhyloNode}
 *
 * @author Michael Scharfe
 *
 */
public class PhyloTree implements Cloneable, Storable {

    private PhyloNode root;
    private String name;

    private final static String XML_TAG = "PHYLO_TREE";
    
    /**
     * Construct an instance of the class PhyloTree
     *
     * @param name the name of the PhyloTree
     * @param root the root node
     */
    public PhyloTree(String name, PhyloNode root) {
    	this.name = name;
    	this.root = root;
    }


    /**
     * The standard constructor for the interface {@link de.jstacs.Storable}.
     * Constructs a {@link PhyloTree} out of an XML representation.
     *
     * @param xml the XML representation as {@link StringBuffer}
     *
     * @throws NonParsableException
     *             if the {@link PhyloTree} could not be reconstructed out of
     *             the {@link StringBuffer} <code>xml</code>
     *
     */
    public PhyloTree(StringBuffer xml) throws NonParsableException {

        xml = XMLParser.extractForTag( xml, "PHYLO_TREE" );
        this.root = (PhyloNode)XMLParser.extractObjectForTags( xml, "root");
        this.name = (String) XMLParser.extractObjectForTags( xml, "name" );
    }

    @Override
    public PhyloTree clone() throws CloneNotSupportedException {
        PhyloTree clone = (PhyloTree) super.clone();
        clone.root = root.clone();
        return clone;
    }

    /**
     * This method returns the name of the PhyloTree
     *
     * @return the name of the tree
     */
    public String getName() {
        return name;
    }

    /**
     * This method returns the root node of the tree
     *
     * @return the root node
     */
    public PhyloNode getRoot() {
        return root;
    }
    
    /**
     * This method returns the total number of nodes in the tree
     *
     * @return the number of nodes in the tree
     */
    public int getNumberOfNodes() {
    	return root.getNumberOfAllNodesBelow() + 1;
    }
    
    /**
     * This method returns a list of {@link PhyloNode}s that represent the leafs of the tree
     *
     * @return the leafs of the tree
     */
    public ArrayList<PhyloNode> getAllLeafs() {
    	return root.getAllLeafs();
    }

    public StringBuffer toXML() {
       StringBuffer xml = new StringBuffer();
       XMLParser.appendObjectWithTags( xml, root, "root" );
       XMLParser.appendObjectWithTags( xml, name, "name" );
       XMLParser.addTags( xml, "PHYLO_TREE" );
       return xml;
    }
}
