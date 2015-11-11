package de.jstacs.clustering.hierachical;

import java.util.Arrays;
import java.util.LinkedList;

import de.jstacs.Storable;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Pair;



public class ClusterTree<T> implements Storable{

	private double distance;
	private T[] elements;
	private ClusterTree<T>[] subTrees;
	private final int originalIndex;
	
	private int[][][] pred;
	
	public ClusterTree(T leaf, int originalIndex){
		this.elements = (T[])ArrayHandler.cast( new Object[]{leaf} );
		this.distance = Double.NEGATIVE_INFINITY;
		this.originalIndex = originalIndex;
	}
	
	public ClusterTree(double distance, int originalIndex, ClusterTree<T>... subTrees){
		this.distance = distance;
		this.subTrees = subTrees.clone();
		this.originalIndex = originalIndex;
		setElements();
	}
	
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
	//returns dmat with new indexes
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
	
	
	public int getOriginalIndex(){
		return originalIndex;
	}
	
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
	
	public ClusterTree<T>[] getSubTrees(){
		return subTrees;
	}
	
	public double getDistance(){
		return distance;
	}
	
	public double getMaximumDistance(){
		return distance;
	}
	
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
	
	public T[] getClusterElements(){
		return elements;
	}

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
	
	public boolean isLeaf(){
		return this.subTrees == null;
	}
	
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
	
	
}
