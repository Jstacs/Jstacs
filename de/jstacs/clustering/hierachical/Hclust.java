package de.jstacs.clustering.hierachical;

import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.clustering.distances.DistanceMetric;
import de.jstacs.io.ArrayHandler;


public class Hclust<T> {

	public enum Linkage{
		SINGLE,
		AVERAGE,
		COMPLETE
	}
	
	private DistanceMetric<T> metric;
	private Linkage linkage;
	
	public Hclust(DistanceMetric<T> metric, Linkage linkage){
		this.metric = metric;
		this.linkage = linkage;
	}
	
	public ClusterTree<T> cluster(T... objects) throws Exception{
		double[][] distMat = DistanceMetric.getPairwiseDistanceMatrix( metric, objects );
		return cluster(distMat, objects);
	}
	
	public ClusterTree<Integer> cluster(double[][] distMat, LinkedList<ClusterTree<Integer>> list, int indexOff){
				
		int oi = -indexOff-1;
		while(list.size() > 1){
			Iterator<ClusterTree<Integer>> it = list.iterator();
			int mini = -1, minj = -1;
			double min = Double.POSITIVE_INFINITY;
			int i=0;
			while(it.hasNext()){
				ClusterTree<Integer> tree = it.next();
				Iterator<ClusterTree<Integer>> it2 = list.iterator();
				for(int j=0;j<i;j++){
					ClusterTree<Integer> tree2 = it2.next();
					double dist = getDistance(distMat, tree,tree2);
					if(dist < min){
						min = dist;
						mini = i;
						minj = j;
					}
				}
				i++;
			}
			
			ClusterTree<Integer> nt = new ClusterTree<Integer>( min, oi, list.get( mini ), list.get( minj ) );
			oi--;
			list.remove( Math.max( mini, minj ) );
			list.remove( Math.min( mini, minj ) );
			list.add( nt );
		}
		
		return list.get( 0 );
	}
	
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
	
	public ClusterTree<T> cluster(double[][] distMat, T... objects){
		LinkedList<ClusterTree<Integer>> list = new LinkedList<ClusterTree<Integer>>();
		for(int i=0;i<objects.length;i++){
			list.add( new ClusterTree<Integer>( i, i ) );
		}
		
		ClusterTree<Integer> tree = cluster(distMat, list, 0);
		
		return createTree( tree, objects );
		
	}

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
	
	public double getDistance( double[][] distMat, ClusterTree<Integer> tree, ClusterTree<Integer> tree2 ) {
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
