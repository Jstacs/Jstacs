/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.clustering.hierachical;

import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.clustering.distances.DistanceMetric;
import de.jstacs.io.ArrayHandler;

/**
 * Class for clustering a set of elements of the same kind (<code>T</code>) hierarchically
 * using agglomerative clustering with different linkage methods.
 * 
 * @author Jan Grau
 *
 * @param <T> the type of the elements clustered
 */
public class Hclust<T> {

	/**
	 * The linkage method for clustering
	 * @author Jan Grau
	 *
	 */
	public enum Linkage{
		/**
		 * Single linkage, i.e., minimum distance between cluster elements
		 */
		SINGLE,
		/**
		 * Average linkage, i.e., average distance between cluster elements (UPGMA)
		 */
		AVERAGE,
		/**
		 * Complete linkage, i.e., maximum distance between cluster elements
		 */
		COMPLETE
	}
	
	private DistanceMetric<T> metric;
	private Linkage linkage;
	
	/**
	 * Creates a new object for clustering using the supplied distance metric and linkage method.
	 * @param metric the distance metric
	 * @param linkage the linkage method
	 */
	public Hclust(DistanceMetric<T> metric, Linkage linkage){
		this.metric = metric;
		this.linkage = linkage;
	}
	
	/**
	 * Clusters the supplied objects and return the resulting cluster tree.
	 * @param objects the objects to be clustered
	 * @return the cluster tree
	 * @throws Exception if the distance matrix could not be created using the {@link DistanceMetric} of this {@link Hclust} object
	 */
	public ClusterTree<T> cluster(T... objects) throws Exception{
		double[][] distMat = DistanceMetric.getPairwiseDistanceMatrix( metric, objects );
		return cluster(distMat, objects);
	}
	
	/**
	 * Further clusters the supplied cluster trees using the given distance matrix
	 * and creating original indexes for the inner node starting at -indexOff-1 in descending order
	 * @param distMat the distance matrix
	 * @param list the list of previous trees
	 * @param indexOff the offset on the original indexes
	 * @return the joint cluster tree
	 */
	public ClusterTree<Integer> cluster(double[][] distMat, LinkedList<ClusterTree<Integer>> list, int indexOff){
				
		int oi = -indexOff-1;
		while(list.size() > 1){
			Iterator<ClusterTree<Integer>> it = list.iterator();
			int mini = -1, minj = -2; //mini>minj
			double min = Double.POSITIVE_INFINITY;
			int i=0;
			while(it.hasNext()){
				ClusterTree<Integer> tree = it.next();
				Iterator<ClusterTree<Integer>> it2 = list.iterator();
				for(int j=0;j<i;j++){
					ClusterTree<Integer> tree2 = it2.next();
					double dist = getDistance(linkage, distMat, tree,tree2);
					if(dist < min){
						min = dist;
						mini = i;
						minj = j;
						//still mini>minj, since i>j
					}
				}
				i++;
			}
			
			ClusterTree<Integer> nt = new ClusterTree<Integer>( min, oi, list.get( mini ), list.get( minj ) );
			oi--;
			list.remove( mini );
			list.remove( minj );
			list.add( nt );
		}
		
		return list.get( 0 );
	}
	
	/**
	 * Clusters the given leaf trees using the supplied distance matrix
	 * @param indexOff the offset on the original indexes
	 * @param distMat the distance matrix
	 * @param leaves the leaves
	 * @return the joint cluster tree
	 * @see Hclust#cluster(int, double[][], ClusterTree[])
	 */
	public ClusterTree<T> cluster(int indexOff, double[][] distMat, ClusterTree<T>[] leaves) {
		LinkedList<ClusterTree<Integer>> list = new LinkedList<ClusterTree<Integer>>();
		T[] objects = (T[])new Object[leaves.length];
		double[][] subMat = new double[leaves.length][leaves.length];
		for(int i=0;i<leaves.length;i++){
			list.add( new ClusterTree<Integer>(i,leaves[i].getOriginalIndex()) );
			objects[i] = leaves[i].getClusterElements()[0];
			for(int j=0;j<leaves.length;j++){
				subMat[i][j] = distMat[leaves[i].getOriginalIndex()][leaves[j].getOriginalIndex()];
			}
		}
		
		ClusterTree<Integer> tree = cluster(subMat,list,indexOff);
		
		return createTree(tree,objects);
	}
	
	/**
	 * Clusters the given objects using the supplied distance matrix, which must be in the same 
	 * order as the elements provided in <code>objects</code>. 
	 * @param distMat the distance matrix
	 * @param objects the objects to the clustered
	 * @return the cluster tree
	 */
	public ClusterTree<T> cluster(double[][] distMat, T... objects){
		LinkedList<ClusterTree<Integer>> list = new LinkedList<ClusterTree<Integer>>();
		for(int i=0;i<objects.length;i++){
			list.add( new ClusterTree<Integer>( i, i ) );
		}
		
		ClusterTree<Integer> tree = cluster(distMat, list, 0);
		
		return createTree( tree, objects );
		
	}

	/**
	 * Cuts the cluster tree at the specified distance and returns the leaf elements
	 * grouped by their origin in the sub-trees below the cut
	 * @param tree the tree
	 * @param distance the cut distance
	 * @param <T> the type of the leaves
	 * @return the leaf elements
	 */
	public static <T> T[][] cutTree(ClusterTree<T> tree, double distance){
		LinkedList<T[]> list = new LinkedList<T[]>();
		fillCutTree( tree, distance, list );
		return (T[][])ArrayHandler.cast( list.toArray( (T[][])new Object[0][] ) );
	}
	
	private static <T> void fillCutTree(ClusterTree<T> tree, double distance, LinkedList<T[]> list){
		if(tree.getDistance() <= distance){
			list.add( tree.getClusterElements() );
		}else{
			ClusterTree<T>[] subs = tree.getSubTrees();
			for(int i=0;i<subs.length;i++){
				fillCutTree( subs[i], distance, list );
			}
		}
	}
	
	/**
	 * Cuts the cluster tree at the given distance and returns the sub-trees below the cut.
	 * @param distance the distance
	 * @param tree the tree
	 * @param <T> the type of the leaves
	 * @return the sub-trees
	 */
	public static <T> ClusterTree<T>[] cutTree(double distance, ClusterTree<T> tree){
		LinkedList<ClusterTree<T>> list = new LinkedList<ClusterTree<T>>();
		fillCutTree( tree, list, distance );
		return list.toArray( (ClusterTree<T>[])new ClusterTree[0] );
	}
	
	private static <T> void fillCutTree(ClusterTree<T> tree, LinkedList<ClusterTree<T>> list, double distance){
		if(tree.getDistance() <= distance){
			list.add( tree );
		}else{
			ClusterTree<T>[] subs = tree.getSubTrees();
			for(int i=0;i<subs.length;i++){
				fillCutTree( subs[i], list, distance );
			}
		}
	}
	
	/**
	 * Creates a cluster tree given an index tree using the original indexes referring to the indexes
	 * of elements in <code>objects</code>.
	 * @param intTree the index tree
	 * @param objects the objects filled into the tree instead of the indexes
	 * @return the cluster tree
	 * @see ClusterTree#getIndexTree()
	 */
	public ClusterTree<T> createTree(ClusterTree<Integer> intTree, T... objects){
		ClusterTree<Integer>[] intSubs = intTree.getSubTrees();
		if(intSubs == null){
			return new ClusterTree<T>( objects[intTree.getClusterElements()[0]], intTree.getOriginalIndex() );
		}else{
			ClusterTree<T>[] subs = new ClusterTree[intSubs.length];
			for(int i=0;i<subs.length;i++){
				subs[i] = createTree( intSubs[i], objects );
			}
			return new ClusterTree<T>( intTree.getDistance(), intTree.getOriginalIndex(), subs );
		}
	}
	
	/**
	 * Returns the distance between the two supplied trees using the linkage method of this {@link Hclust} object
	 * and the given distance matrix.
	 * @param distMat the distance matrix
	 * @param tree the first tree
	 * @param tree2 the second tree
	 * @return the distance between the trees
	 */
	public static double getDistance( Linkage linkage, double[][] distMat, ClusterTree<Integer> tree, ClusterTree<Integer> tree2 ) {
		double dist = linkage == Linkage.SINGLE ? Double.POSITIVE_INFINITY : ( linkage == Linkage.AVERAGE ? 0 : Double.NEGATIVE_INFINITY );
		
		Integer[] el1 = tree.getClusterElements();
		Integer[] el2 = tree2.getClusterElements();
		for(int i=0;i<el1.length;i++){
			for(int j=0;j<el2.length;j++){
				double d = distMat[ Math.max( el1[i], el2[j] ) ][ Math.min( el1[i], el2[j] ) ];
				if(linkage == Linkage.SINGLE){
					if(d < dist){
						dist = d;
					}
				}else if(linkage == Linkage.AVERAGE){
					dist += d / (el1.length*el2.length);
				}else if(linkage == Linkage.COMPLETE){
					if(d > dist){
						dist = d;
					}
				}else{
					throw new RuntimeException( "Linkage not supported" );
				}
			}
		}
		return dist;
	}
	
}
