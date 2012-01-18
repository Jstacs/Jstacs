package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.MultiDimensionalDiscreteSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.SamplingEmission;
import de.jstacs.sequenceScores.statisticalModels.trainable.phylo.PhyloNode;
import de.jstacs.sequenceScores.statisticalModels.trainable.phylo.PhyloTree;
import de.jstacs.utils.Normalisation;
import java.util.*;
import javax.naming.OperationNotSupportedException;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import java.util.ArrayList;

/**
 * Phylogenetic discrete emission
 * This class uses a phylogenetic tree to describe multidimensional data
 * It implements Felsensteins model for nucleotide substitution (F81)
 *
 * @author Michael Scharfe
 */
public class PhyloDiscreteEmission extends DiscreteEmission implements SamplingEmission{

    private PhyloTree tree;
    private double[][] pruningMatrix;
    private HashMap<Integer, ArrayList<PhyloNode>> treeLayers;
    private double logProbFromStatistic;

    /**
     * This is a simple constructor for a {@link PhyloDiscreteEmission} based on the equivalent sample size.
     *
     * @param con the {@link AlphabetContainer} of this emission
     * @param ess the equivalent sample size (ess) of this emission that is equally distributed over all parameters
     * @param t the phylogenetic tree {@link PhyloTree}
     */
    public PhyloDiscreteEmission( AlphabetContainer con, double ess , PhyloTree t) {
            this( con, getHyperParams( ess, (int) con.getAlphabetLengthAt( 0 ) ), t );
    }

    /**
     * 
     * This is a simple constructor for a {@link DiscreteEmission} defining the individual hyper parameters.
     *
     * @param con the {@link AlphabetContainer} of this emission
     * @param hyperParams the individual hyper parameters for each parameter
     * @param t the phylogenetic tree {@link PhyloTree}
     */
    public PhyloDiscreteEmission( AlphabetContainer con, double[] hyperParams, PhyloTree t ) {
        super(con,hyperParams);
        this.tree = t;
        setTreeLayers();
    }

   
    @Override
    public PhyloDiscreteEmission clone() throws CloneNotSupportedException {
        PhyloDiscreteEmission clone = (PhyloDiscreteEmission) super.clone();
        clone.tree = tree.clone();
        if(pruningMatrix != null){
            clone.pruningMatrix = ArrayHandler.clone( pruningMatrix );
        }
        clone.logProbFromStatistic = logProbFromStatistic;
        clone.treeLayers = (HashMap<Integer, ArrayList<PhyloNode>>) treeLayers.clone();
        return clone;
    }

    private static double[] getHyperParams( double ess, int number ){
        double[] res = new double[number];
        Arrays.fill( res, ess/(double) number );
        return res;
    }

    @Override
    public double getLogProbFor(boolean forward, int startPos, int endPos, Sequence seq) throws OperationNotSupportedException {

        //TODO forward==F (currently only forward==T)

        for(int pos = startPos; pos <= endPos; pos++) {
            fillPruningMatrix(pos, seq);
        }   

        return Normalisation.getLogSum(pruningMatrix[0]);
    }



    @Override
    public double getLogProbAndPartialDerivationFor( boolean forward, int startPos, int endPos,
			IntList indices, DoubleList partDer, Sequence seq) throws OperationNotSupportedException {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void addGradientOfLogPriorTerm(double[] gradient, int offset) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void estimateFromStatistic() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    private static final String XML_TAG = "PhyloDiscreteEmission";


    @Override
    protected void appendFurtherInformation(StringBuffer xml) {
        XMLParser.appendObjectWithTags( xml, pruningMatrix, "pruningMatrix" );
        XMLParser.appendObjectWithTags( xml, logProbFromStatistic, "logProbFromStatistic" );
        XMLParser.appendObjectWithTags( xml, tree, "PhyloTree" );
    }

    @Override
    protected void extractFurtherInformation(StringBuffer xml) throws NonParsableException {
        pruningMatrix = (double[][]) XMLParser.extractObjectForTags( xml, "pruningMatrix" );
        logProbFromStatistic = (Double) XMLParser.extractObjectForTags( xml, "logProbFromStatistic" );
        tree = (PhyloTree) XMLParser.extractObjectForTags( xml, "PhyloTree" );
        setTreeLayers();
    }

    
    /**
     * This methods fills the pruning-matrix for efficient likelihood calculation
     * It implements the Felsenstein model for dna-substitution (F81)
     *
     * @param pos The position in the sequence
     * @param seq The (multidimensional) sequence
     */
    private void fillPruningMatrix(int pos, Sequence seq) {

        HashMap<Integer, PhyloNode> hm = new HashMap<Integer, PhyloNode>();
        int alphabetSize = (int)(con.getMaximalAlphabetLength());
        pruningMatrix = new double[tree.getNumberOfNodes()][alphabetSize];
        int[] alignment = (int[])(seq.getEmptyContainer());
        seq.fillContainer(alignment, pos);
        SequenceAnnotation[][] annot = ((MultiDimensionalDiscreteSequence)seq).getAnnotations();
        ArrayList<PhyloNode> tmpNodes, children;
        double[] tmp = new double[alphabetSize];

        //bottom -> up
        for(int l = treeLayers.size()-1; l>=0; l--) {
            tmpNodes = treeLayers.get(l);

            for(PhyloNode n : tmpNodes) {
                children = n.getChildrenNodes();

                //inner nodes
                if(children.size() > 0) {
                    for(int a = 0; a < alphabetSize; a++) {
                        for(int b = 0; b < alphabetSize; b++) {
                            tmp[b] = (a==b) ? Math.log(1d-n.getWeight() + n.getWeight()*probs[0][a]) : Math.log(n.getWeight()*probs[0][a]);
                            
                            for(int c = 0; c < children.size(); c++) {
                                tmp[b] += pruningMatrix[children.get(c).getId()][b];
                            }
	                }
	                pruningMatrix[n.getId()][a] += Normalisation.getLogSum(tmp);
                    }
                }
                //leafs
                else {
                    for(int s = 0; s < alignment.length; s++) {
                        if(annot[s][0].getIdentifier().equals(n.getName())) {
                            for(int a = 0; a < alphabetSize; a++) {
                                pruningMatrix[n.getId()][a] = alignment[s] == a ? Math.log(1d-n.getWeight() + n.getWeight()*probs[0][a]) : Math.log(n.getWeight()*probs[0][a]);
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Returns the log posterior of the proposal distribution for the current statistic
     *
     * @return the log posterior
     */
    public double getLogProposalPosteriorFromStatistic() {
        return super.getLogPosteriorFromStatistic();
    }

    @Override
    public double getLogPosteriorFromStatistic() {
        
        return logProbFromStatistic + getLogPriorTerm();
    }

    @Override
    public void resetStatistic() {
        super.resetStatistic();
        logProbFromStatistic = 0d;
    }



    @Override
    public void addToStatistic(boolean forward, int startPos, int endPos, double weight, Sequence seq) throws OperationNotSupportedException {
        int s, e;
        Sequence current;
        int[] alignment = (int[])(seq.getEmptyContainer());
        
        if( forward ) {
                current = seq;
                s = startPos;
                e = endPos;
        } else {
                current = seq.reverseComplement();
                int len = current.getLength();
                s = len - endPos -1;
                e = len - startPos -1;
        }
              
        while( s <= e ) {
            current.fillContainer(alignment, s);
            for(int i = 0; i < alignment.length; i++) {
                statistic[0][alignment[i]] += weight;
            }
            s++;
        }
        logProbFromStatistic += getLogProbFor(forward, startPos, endPos, seq);
    }  

    /* Maps phylogenetic tree to HashMap
     *
     */
    private void setTreeLayers() {

        treeLayers = new HashMap<Integer, ArrayList<PhyloNode>>();
        ArrayList<PhyloNode> nodeList = new ArrayList<PhyloNode>();
        PhyloNode root = tree.getRoot();
        root.setWeight(1d);
        nodeList.add(root);
        int layer = 0;

        while(nodeList.size() > 0) {
            treeLayers.put(layer, nodeList);
            nodeList = new ArrayList<PhyloNode>();
            for(PhyloNode n : treeLayers.get(layer)) {
                nodeList.addAll(n.getChildrenNodes());
            }
            layer++;
        }
    }

    @Override
    public double getLogGammaScoreFromStatistic() {
       throw new UnsupportedOperationException("Not supported yet.");
    }

}
