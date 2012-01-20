package de.jstacs.sequenceScores.statisticalModels.trainable.phylo;

import java.util.ArrayList;
import java.util.Arrays;

import de.jstacs.Storable;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * This class implements a node in a {@link PhyloTree}
 * A PhyloNode contains some basic informations of itself and the incoming edge
 * Furthermore it contains a list of {@link PhyloNode}s that represent the children nodes
 *
 * @author Michael Scharfe
 *
 */
public class PhyloNode implements Cloneable, Storable {

    private String name;
    private double weight;
    private ArrayList<PhyloNode> children;
    private int id;

    private final static String XML_TAG = "PHYLO_NODE";

    /**
     * Basic constructor
     *
     */
    public PhyloNode() {
        children = new ArrayList<PhyloNode>();
        id = Integer.MIN_VALUE;
        weight = Double.MIN_VALUE;
    }

    /**
     * The standard constructor for the interface {@link de.jstacs.Storable}.
     * Constructs a {@link PhyloNode} out of an XML representation.
     *
     * @param xml the XML representation as {@link StringBuffer}
     *
     * @throws NonParsableException
     *             if the {@link PhyloNode} could not be reconstructed out of
     *             the {@link StringBuffer} <code>xml</code>
     *
     */
    public PhyloNode(StringBuffer xml) throws NonParsableException {
        
        xml = XMLParser.extractForTag( xml, "PHYLO_NODE" );
        this.name = (String) XMLParser.extractObjectForTags( xml, "name" );
        this.weight = (Double) XMLParser.extractObjectForTags( xml, "weight" );

        children = new ArrayList<PhyloNode>();

        if(XMLParser.extractObjectForTags( xml, "numberOfChildren", Integer.class) > 0) {
            PhyloNode[] c = XMLParser.extractObjectForTags( xml, "children", PhyloNode[].class);
            if(c.length>0){
            	this.children.addAll(Arrays.asList(c));
            }
        }
        this.id = (Integer) XMLParser.extractObjectForTags( xml, "id" );
    }

    @Override
    public PhyloNode clone() throws CloneNotSupportedException {
        PhyloNode clone = (PhyloNode) super.clone();
        ArrayList<PhyloNode> c = new ArrayList<PhyloNode>(children.size());
        for(PhyloNode item: children) c.add(item.clone());
        clone.children = c;
        return clone;
    }

    /**
     * This method adds a children node to the current instance
     *
     * @param node the children node
     */
    public void addChild(PhyloNode node) {
    	children.add(node);
    }
    
    /**
     * This method set a name for the current instance
     *
     * @param name the name of the PhyloNode
     */
    public void setName(String name) {
    	this.name = name;
    }

    /**
     * This method returns the name of the current instance
     *
     * @return the name of the PhyloNode
     */
    public String getName() {
        return name;
    }

    /**
     * This method set the ID of the current PhyloNode
     * The ID should be unique in the PhyloTree
     *
     * @param id the ID of the node
     */
    public void setId(int id) {
        this.id = id;
    }

    /**
     * This method returns the ID of the current PhyloNode
     * 
     * @return the ID of the node
     */
    public int getId() {
        return id;
    }

    /**
     * This method set the weight (length, rate ...) for the incoming edge
     *
     * @param w the weight
     */
    public void setWeight(double w) {
    	weight = w;
    }
    
    /**
     * This method return the weight (length, rate ...) for the incoming edge
     *
     * @return the weight
     */
    public double getWeight() {
    	return weight;
    }
    
    /**
     * This method returns the total number of {@link PhyloNode}s in the subtree starting from this instance
     *
     * @return the total number of nodes
     */
    public int getNumberOfAllNodesBelow() {
    	
    	int n = 0;
    	
    	for(PhyloNode node : children) {
    		n += node.getNumberOfAllNodesBelow() + 1;
    	}
    	return n;
    }

    /**
     * This method returns a list of {@link PhyloNode}s that are children of this instance
     *
     * @return the children nodes
     */
    public ArrayList<PhyloNode> getChildrenNodes() {
        return children;
    }
    
    /**
     * This method returns a list of {@link PhyloNode}s that are leafs in the subtree starting from this instance
     *
     * @return the list of leafs
     */
    public ArrayList<PhyloNode> getAllLeafs() {
    	
    	ArrayList<PhyloNode> l = new ArrayList<PhyloNode>();
    	
    	if(children.isEmpty())
    		l.add(this);
    	
    	for(PhyloNode n : children) {
    		l.addAll(n.getAllLeafs());
    	}
    	
    	return l;
    }
    
    public StringBuffer toXML() {
       StringBuffer xml = new StringBuffer();
       XMLParser.appendObjectWithTags( xml, name, "name" );
       XMLParser.appendObjectWithTags( xml, weight, "weight" );

       XMLParser.appendObjectWithTags( xml, children.size(), "numberOfChildren");
       if(children.size() > 0) {
           PhyloNode[] c = children.toArray( new PhyloNode[0] );
           XMLParser.appendObjectWithTags( xml, c, "children");
       }
      
       XMLParser.appendObjectWithTags( xml, id, "id" );
       XMLParser.addTags( xml, "PHYLO_NODE" );
       return xml;
    }
    
    
}

