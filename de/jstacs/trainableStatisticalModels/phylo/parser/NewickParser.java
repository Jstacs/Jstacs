package de.jstacs.trainableStatisticalModels.phylo.parser;

import de.jstacs.NonParsableException;
import de.jstacs.trainableStatisticalModels.phylo.PhyloNode;
import de.jstacs.trainableStatisticalModels.phylo.PhyloTree;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.util.EmptyStackException;
import java.util.Stack;



/**
 * This class implements a simple newick parser and allows the construction of a
 * {@link PhyloTree}
 *
 * @author Michael Scharfe
 */
public class NewickParser {

    private StreamTokenizer tokenizer;
    private PhyloNode rootNode;
    
    /**
     * This constructor initialize the newick parser
     *
     * @param b the {@link BufferedReader} which reads the newick file
     *
     */
    public NewickParser(BufferedReader b)
    {
        tokenizer = new StreamTokenizer(b);
        tokenizer.eolIsSignificant(false);
        tokenizer.quoteChar('"');
        tokenizer.wordChars('\'', '\''); // quote problem, turn this into a prime symbol?
        tokenizer.wordChars('!', '!'); // 33
        tokenizer.wordChars('#', '&'); // 35-38
        tokenizer.wordChars('*', '+'); // 42-43
        tokenizer.wordChars('-', '/'); // 45-47
        tokenizer.wordChars('<', '<'); // 60
        tokenizer.wordChars('>', '@'); // 62-64
        tokenizer.wordChars('^', '`'); // 93-96
        tokenizer.wordChars('{', '~'); // 123-126
    }
	
        /**
         * This method construct a {@link PhyloTree} from the given input stream
         *           
         * @return the {@link PhyloTree}
         * @throws NonParsableException if the stream is not parsable
         */
        public PhyloTree tokenize() throws NonParsableException {

        final char openBracket = '(', 
        	closeBracket = ')', 
        	childSeparator = ',',
        	infoSeparator = ':';
        int thisToken;
        boolean EOT = false;
        boolean nameNext = true;
        PhyloNode lastNamed = null, newNode;
        Stack nodeStack = new Stack();
        int id = 0;

        rootNode = new PhyloNode();
        rootNode.setId(id++);
        PhyloTree tree = new PhyloTree(null,rootNode);	//TODO parse tree name
        nodeStack.push(rootNode);
        
        try
        {
            while (EOT == false &&
                    (thisToken = tokenizer.nextToken()) != StreamTokenizer.TT_EOF)
            {
	            switch (thisToken)
	            {
	            	case StreamTokenizer.TT_WORD:
	            	    if (!nameNext)
	            	        System.err.println("Error: didn't expect this name here: " + tokenizer.sval);
	            	    lastNamed = popNodeAndName(tokenizer.sval, nodeStack);
	            		nameNext = false;
	            		break;
	            	case StreamTokenizer.TT_NUMBER:
	            		if (nameNext)	//new node name
	            		    lastNamed = popNodeAndName(tokenizer.sval, nodeStack);
	            		else			 
	            		{				//weight of last node
	            		    if (lastNamed != null)
	            		        lastNamed.setWeight(tokenizer.nval);
	            		    else
	            		        System.err.println("Error: can't set value " + tokenizer.nval + " to a null node");
	            		    lastNamed = null;
	            		}
	            		nameNext = false;
	            		break;
	            	case infoSeparator:
	            	    if (nameNext)	//last node has no name
	            	        lastNamed = popNodeAndName(null, nodeStack);
	            	    nameNext = false;
	            	    break;
	            	case StreamTokenizer.TT_EOF:
	            	    if (nameNext)
	            	        lastNamed = popNodeAndName(null, nodeStack);
	            	    EOT = true;
	            	    nameNext = false;
	            	    break;
	            	case openBracket:	//new node
                            newNode = new PhyloNode();
                            newNode.setId(id++);
	            	    nodeStack.push(newNode);
	            	    nameNext = true;
	            	    break;
	            	case closeBracket:
	            	    if (nameNext)
	            	        lastNamed = popNodeAndName(null, nodeStack);
	            	    nameNext = true;
	            	    break;
	            	case childSeparator:
	            	    if (nameNext)
	            	        lastNamed = popNodeAndName(null, nodeStack);
                            newNode = new PhyloNode();
                            newNode.setId(id++);
	            	    nodeStack.push(newNode);
	            	    nameNext = true;
	            	    break;
	            	default:
	            	    //TODO
	            		break;
	            }
	        }
	    }
        catch (IOException e)
        {
            NonParsableException npe = new NonParsableException( e.getMessage() );
            throw npe;
        }
	return tree;
    }
	
	
    private PhyloNode popNodeAndName(String name, Stack nodeStack) throws NonParsableException
    {
	    PhyloNode node = (PhyloNode)nodeStack.pop();
	    node.setName(name==null ? "" : name);

	    try
	    {
	    	PhyloNode parent = (PhyloNode) nodeStack.peek();
	    	parent.addChild(node);
	    }
	    catch (EmptyStackException e)
	    {
                NonParsableException npe = new NonParsableException( e.getMessage() );
                throw npe;
	    }
	    return node;
    }
}
