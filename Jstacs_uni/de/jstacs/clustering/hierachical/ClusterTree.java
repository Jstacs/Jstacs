package de.jstacs.clustering.hierachical;

import java.util.Arrays;
import java.util.LinkedList;

import de.jstacs.Storable;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Pair;


/**
 * Class for a generic cluster tree with leaves of type <code>T</code>.
 * Cluster trees for a given set of leaf elements may be obtained using {@link Hclust}.
 * 
 * @author Jan Grau
 *
 * @param <T> the type of the leaves
 */
public class ClusterTree<T> implements Storable{

	private double distance;
	private T[] elements;
	private ClusterTree<T>[] subTrees;
	private int originalIndex;
	
	private int[][][] pred;
	
	/**
	 * Creates a new cluster tree for a given leaf element (i.e., the tree comprises just this leaf) with the
	 * supplied index in the set of cluster elements
	 * @param leaf the leaf element
	 * @param originalIndex the index of the leaf element in the complete set of cluster elements
	 */
	public ClusterTree(T leaf, int originalIndex){
		this.elements = (T[])ArrayHandler.cast( new Object[]{leaf} );
		this.distance = Double.NEGATIVE_INFINITY;
		this.originalIndex = originalIndex;
	}
	
	/**
	 * Creates a new cluster tree with supplied sub-trees and given distance.
	 * The original index may be a virtual index that helps to identify this inner (or root) node
	 * of the cluster tree later on.
	 * @param distance the distance between the sub-trees
	 * @param originalIndex the original index of this inner (or root) node
	 * @param subTrees the sub-trees
	 */
	public ClusterTree(double distance, int originalIndex, ClusterTree<T>... subTrees){
		this.distance = distance;
		this.subTrees = subTrees.clone();
		this.originalIndex = originalIndex;
		setElements();
	}
	
	/**
	 * Creates a  cluster tree from its XML representation
	 * @param xml the XML representation
	 * @throws NonParsableException if XML could not be parsed
	 */
	public ClusterTree(StringBuffer xml) throws NonParsableException{
		xml = XMLParser.extractForTag( xml, "ClusterTree" );
		originalIndex = (Integer)XMLParser.extractObjectForTags( xml, "originalIndex", Integer.class );
		distance = (Double)XMLParser.extractObjectForTags( xml, "distance" );
		subTrees = (ClusterTree<T>[])XMLParser.extractObjectForTags( xml, "subTrees" );
		if(subTrees == null){
			elements = (T[])XMLParser.extractObjectForTags( xml, "elements" );
		}else{
			setElements();
		}
		
	}
	
	/**
	 * Returns a cluster tree with identical structure as this cluster tree but with all leaves replaced by
	 * integer leaves holding the corresponding original indices.
	 * @return the index cluster tree
	 * @see ClusterTree#ClusterTree(Object, int)
	 */
	public ClusterTree<Integer> getIndexTree(){
		if(subTrees == null){
			return new ClusterTree<Integer>(originalIndex,originalIndex);
		}else{
			ClusterTree<Integer>[] newSubs = new ClusterTree[subTrees.length];
			for(int i=0;i<subTrees.length;i++){
				newSubs[i] = subTrees[i].getIndexTree();
			}
			return new ClusterTree<Integer>(distance,originalIndex,newSubs);
		}
	}
	
	//algo from Ziv Bar-Joseph et al. K-ary Clustering with Optimal Leaf Ordering for Gene Expression Data.
	//dmat: original dmat when building the tree (original index!)
	/**
	 * Orders the leaves of this cluster tree such that adjacent nodes have minimal distance.
	 * Implements the algorithm of Ziv Bar-Joseph et al. "K-ary Clustering with Optimal Leaf Ordering for Gene Expression Data"
	 * 
	 * @param dmat the distance matrix that has also been used to build the tree (using the same indexes)
	 */
	public void leafOrder(double[][] dmat){
		
		Pair<double[][],int[]> pair = forward( dmat );
		double[][] M = pair.getFirstElement();
		
		double min = Double.POSITIVE_INFINITY;
		int mini = -1, minj = -1;
		for(int i=0;i<M.length;i++){
			for(int j=0;j<M[i].length;j++){
				if(M[i][j] < min){
					min = M[i][j];
					mini = i;
					minj = j;
				}
			}
		}
		
		
		backtrace(mini,minj);
		
		/*ClusterTree<T>[] leaves = getLeaves();
		
		double[][] dmat2 = new double[dmat.length][dmat.length];
		for(int i=0;i<leaves.length;i++){
			for(int j=0;j<leaves.length;j++){
				dmat2[i][j] = dmat[ leaves[i].getOriginalIndex() ][ leaves[j].getOriginalIndex() ];
			}
		}
		for(int i=0;i<leaves.length;i++){
			leaves[i].originalIndex = i;
		}
		
		return dmat2;*/
		
		//return dmat;
		
	}
	
	private boolean backtrace(int mini, int minj){
		
		//we are at a leaf
		if(pred == null){
			return false;
		}
		
		//predecessor of left (after turn (!) if necessary)
		int predLeft = pred[0][mini][minj];
		// and right child
		int predRight = pred[1][mini][minj];
		
		//if we have been at the lower triangular matrix of M
		//we need to turn the current tree, i.e., exchange the two children
		boolean turn = mini > minj; 
		
		//offset on larger index
		int oneoff = this.subTrees[0].getNumberOfElements();
		
		if(turn){
			
			//correct index (mini has been larger one)
			mini -= oneoff;
			
			//exchange indexes, because order after turn, and we 
			//haven't turned, yet
			int temp = mini;
			mini = predRight;
			predRight = temp;
			
			temp = minj;
			minj = predLeft;
			predLeft = temp;
		}else{
			//no turn, only correct minj
			minj -= oneoff;
		}
		
		//indexes in subtrees (for M): (mini,predLeft), (predRight,minj)
		boolean turnedLeft = this.subTrees[0].backtrace( mini, predLeft );
		boolean turnedRight = this.subTrees[1].backtrace( predRight, minj );
		
		//turn, if necessary
		if(turn){
			this.reverseOrder();
		}else if(turnedLeft || turnedRight){
			//if children have been turned,
			//we need to update the order of our children
			this.setElements();
		}
		
		//remove to free memory
		this.pred = null;
		
		//return whether any (sub-) tree below this one has been turned
		return turn || turnedLeft || turnedRight;
		
	}
	
	private Pair<double[][],int[]> forward(double[][] dmat){
		//we are a leaf
		if(this.subTrees == null){
			pred = null;
			return new Pair<double[][],int[]>(new double[][]{{0}},new int[]{this.originalIndex});
		}else{
			Pair<double[][],int[]> pairLeft = this.subTrees[0].forward( dmat );
			Pair<double[][],int[]> pairRight = this.subTrees[1].forward( dmat );
			//w and x are the equivalents to our M (DP matrix) for nodes under 
			// the left and right child
			double[][] w = pairLeft.getFirstElement();
			double[][] x = pairRight.getFirstElement();
			//we also need the original indexes of the leaves
			//for getting the right distances from dmat
			int[] leftIdx = pairLeft.getSecondElement();
			int[] rightIdx = pairRight.getSecondElement();
			
			
			//our new matrix of optimal values for given bordering leaves (see paper)
			//matrix will be rather sparse, could be optimized
			double[][] M = new double[this.getNumberOfElements()][this.getNumberOfElements()];
			//pred is the field for backtracking, one for the left and one for the right child
			//(first index), index values in pred refer to indexes in w and x (but not the original indexes)
			pred = new int[2][M.length][M.length];
			for(int i=0;i<M.length;i++){
				Arrays.fill( M[i], Double.POSITIVE_INFINITY );
				Arrays.fill( pred[0][i], -1 );
				Arrays.fill( pred[1][i], -1 );
			}
			
			
			// we collect the original indexes from our children
			int[] origIdx = new int[M.length];
			System.arraycopy( leftIdx, 0, origIdx, 0, leftIdx.length );
			System.arraycopy( rightIdx, 0, origIdx, leftIdx.length, rightIdx.length );
			
			//offset on all indexes from right child, necessary to have a common (but sparse)
			//matrix for all children, altough only combinations of leaves from distinct
			//children are allowed
			int oneoff = this.subTrees[0].getNumberOfElements();
			
			//fill DP matrix (see paper)
			for(int i=0;i<w.length;i++){
				//first minimization
				double[] tempil = new double[x.length];
				int[] minil = new int[x.length];
				for(int l=0;l<tempil.length;l++){
					tempil[l] = Double.POSITIVE_INFINITY;
					
					for(int h=0;h<w[i].length;h++){
						double temp = w[i][h] + dmat[ leftIdx[h] ][ rightIdx[l] ];
						//here, we minimize
						if(temp < tempil[l]){
							tempil[l] = temp;
							minil[l] = h;
						}
					}
				}
				//second minimization
				for(int l=0;l<x.length;l++){
					for(int j=0;j<x[l].length;j++){
						double temp = tempil[l] + x[l][j];
						if(temp < M[i][oneoff + j]){
							M[ i ][ oneoff + j ] = temp;
							M[ oneoff + j ][ i ] = temp;
							
							//indexes for backtracking
							pred[0][ i ][ oneoff + j ] = minil[l];
							pred[1][ i ][ oneoff + j ] = l;
							//if swapped, we also swap indexes
							pred[0][ oneoff + j ][ i ] = l;
							pred[1][ oneoff + j ][ i ] = minil[l];
						}
					}
				}
			}
			
			/*System.out.println( "Leaves under "+this.originalIndex );
			System.out.println(Arrays.toString( this.getLeaves() ) );
			System.out.println("Original indexes:");
			System.out.println(Arrays.toString( origIdx ));
			
			System.out.println("M for "+this.originalIndex);
			for(int i=0;i<M.length;i++){
				for(int j=0;j<M[i].length;j++){
					System.out.print(M[i][j]+"\t");
				}
				System.out.println();
			}
			System.out.println("Pred for "+this.originalIndex);
			for(int i=0;i<pred.length;i++){
				for(int j=0;j<pred[i].length;j++){
					for(int k=0;k<pred[i][j].length;k++){
						System.out.print(pred[i][j][k]+"\t");
					}
					System.out.println();
				}
				System.out.println();
			}
			*/
			
			return new Pair<double[][],int[]>(M,origIdx);
			
		}
	}
	
	/**
	 * Returns the original index of the root node of this cluster tree
	 * @return the original index
	 */
	public int getOriginalIndex(){
		return originalIndex;
	}
	
	/**
	 * Reverses the order of the child trees of this cluster tree root node.
	 */
	public void reverseOrder(){
		//System.out.println("reversing "+originalIndex);
		if(subTrees == null || subTrees.length == 1){
			return;
		}
		ClusterTree<T>[] temp = new ClusterTree[subTrees.length];
		for(int i=0;i<subTrees.length;i++){
			temp[i] = subTrees[subTrees.length-1-i];
		}
		subTrees = temp;
		//System.out.println("before: "+Arrays.toString( elements ));
		setElements();
		//System.out.println("after: "+Arrays.toString( elements ));
	}
	
	@Override
	public StringBuffer toXML() {
		StringBuffer sb = new StringBuffer();
		XMLParser.appendObjectWithTags( sb, originalIndex, "originalIndex" );
		XMLParser.appendObjectWithTags( sb, distance, "distance" );
		XMLParser.appendObjectWithTags( sb, subTrees, "subTrees" );
		if(subTrees == null){
			XMLParser.appendObjectWithTags( sb, elements, "elements" );
		}
		XMLParser.addTags( sb, "ClusterTree" );
		return sb;
	}
	
	/**
	 * Sets the copy references of the leave nodes of this cluster
	 * tree to the elements of its leaves in the current order.
	 */
	public void setElements(){
		int num = 0;
		for(int i=0;i<subTrees.length;i++){
			num += subTrees[i].elements.length;
		}
		Object[] elements = new Object[num];
		for(int i=0,k=0;i<subTrees.length;i++){
			for(int j=0;j<subTrees[i].elements.length;j++,k++){
				elements[k] = subTrees[i].elements[j];
			}
		}
		this.elements = (T[]) ArrayHandler.cast( elements );
	}
	
	/**
	 * Returns the sub-trees of this cluster tree root node
	 * @return the sub-trees
	 */
	public ClusterTree<T>[] getSubTrees(){
		return subTrees;
	}
	
	/**
	 * Returns the distance between the child trees of this cluster tree root node.
	 * @return the distance
	 */
	public double getDistance(){
		return distance;
	}
	
	/**
	 * Returns the maximum distance of trees under this root node.
	 * In this implementation identical to {@link ClusterTree#getDistance()}.
	 * Basically used for plotting the tree structure.
	 * @return the distance
	 */
	public double getMaximumDistance(){
		return distance;
	}
	
	/**
	 * Returns the minimum distance of trees under this root node.
	 * Basically used for plotting the tree structure.
	 * @return the distance
	 */
	public double getMinimumDistance(){
		if(getNumberOfElements() == 1){
			return 0;
		}else{
			double min = this.getDistance();
			for(int i=0;i<subTrees.length;i++){
				double temp = subTrees[i].getMinimumDistance();
				if(temp < min){
					min = temp;
				}
			}
			return min;
		}
	}
	
	/**
	 * Returns the elements at all leaves in this cluster tree, in the order
	 * of the leaves, from left to right.
	 * @return the elements
	 */
	public T[] getClusterElements(){
		return elements;
	}

	/**
	 * Returns the number of leaves in this cluster tree.
	 * @return the number of leaves
	 */
	public int getNumberOfElements(){
		return elements.length;
	}
	
	public String toString(){
		if(subTrees == null){
			return Arrays.toString( elements );
		}else{
			StringBuffer sb = new StringBuffer();
			sb.append("dist: "+distance+"\n");
			for(int i=0;i<subTrees.length;i++){
				sb.append(i+"( ["+originalIndex+"]\n");
				sb.append(subTrees[i].toString()+"\n");
				sb.append(i+")\n");
			}
			return sb.toString();
		}
	}

	/**
	 * Returns the minimum original index in this cluster tree.
	 * @return the minimum original index
	 */
	public int getMinimumOriginalIndex(){
		if(this.subTrees == null){
			return this.originalIndex;
		}else{
			int min = this.originalIndex;
			for(int i=0;i<subTrees.length;i++){
				if(subTrees[i].getOriginalIndex() < min){
					min = subTrees[i].getOriginalIndex();
				}
			}
			return min;
		}
	}
	
	/**
	 * Returns <code>true</code> if this cluster tree comprises just a leaf.
	 * @return if this tree is a leaf
	 */
	public boolean isLeaf(){
		return this.subTrees == null;
	}
	
	/**
	 * Returns all leaves of this cluster tree as {@link ClusterTree} objects comprising just the corresponding
	 * leaf element
	 * @return the leaf trees
	 */
	public ClusterTree<T>[] getLeaves() {
		if(this.subTrees == null){
			return new ClusterTree[]{this};
		}else{
			LinkedList<ClusterTree<T>> list = new LinkedList<ClusterTree<T>>();
			for(int i=0;i<subTrees.length;i++){
				ClusterTree<T>[] temp = subTrees[i].getLeaves();
				for(int j=0;j<temp.length;j++){
					list.add( temp[j] );
				}
			}
			return list.toArray( new ClusterTree[0] );
		}
	}

	/**
	 * Removes all sub-trees below the inner nodes identified by the original indexes supplied
	 * and creates new leaf nodes instead, which obtain the supplied leaf elements. The modified tree
	 * is returned.
	 * @param rootOriginalIndexes the original indexes of the inner nodes to be replaced
	 * @param newElements the leaf element replacements
	 * @return the modified tree
	 */
	public <S> ClusterTree<S> dropBelow( IntList rootOriginalIndexes, S[] newElements ) {
		int idx = rootOriginalIndexes.contains( originalIndex ); 
		if(idx > -1){
			S el = newElements[idx];
			return new ClusterTree<S>(el,originalIndex);
		}else{
			ClusterTree<S>[] newSubs = new ClusterTree[subTrees.length];
			for(int i=0;i<subTrees.length;i++){
				newSubs[i] = subTrees[i].dropBelow(rootOriginalIndexes, newElements);
			}
			return new ClusterTree<S>(this.getDistance(),originalIndex,newSubs);
		}
	}

	/**
	 * Returns a string representation of this cluster tree in a pseudo newick format.
	 * @return the string representation
	 */
	public String toNewick() {
		return this.toNewick("");
	}

	private String toNewick(String indent) {
		StringBuffer sb = new StringBuffer();
		if(this.subTrees == null){
			sb.append(indent+"("+elements[0].toString()+")\n");
		}else{
			sb.append(indent+this.distance+" (\n");
			for(int i=0;i<subTrees.length;i++){
				sb.append(subTrees[i].toNewick(indent+"\t")+"\n");
			}
			sb.append(indent+")\n");
		}
		return sb.toString();
	}

	/**
	 * Sets the original index (e.g., if elements have been removed from the tree) referring to indexes
	 * in the distance matrix that has been used to build a tree.
	 * @param used the new index
	 */
	public void setOriginalIndex(int used) {
		this.originalIndex = used;
	}
	
	
}
