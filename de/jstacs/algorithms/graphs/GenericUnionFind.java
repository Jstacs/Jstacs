package de.jstacs.algorithms.graphs;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * 
 * @author Jens Keilwagen
 */
public class GenericUnionFind<T> {

	ArrayList<T> index;
	HashMap<T,int[]> uf;
	
	public GenericUnionFind() {
		uf = new HashMap<T, int[]>();
		index = new ArrayList<T>();
	}
	
	public void union( T first, T second ) {
		int firstRoot = find(first);
		int secondRoot = find(second);
		if( firstRoot != secondRoot ) {
			T help=index.get(firstRoot);
			int[] value1 = uf.get(help);
			help=index.get(secondRoot);
			int[] value2 = uf.get(help);

			if( value1[1] > value2[1] ) {
				int[] xxx = value1;
				value1=value2;
				value2=xxx;
			}
			//value1 is the larger component (value1[1]<value2[1])
			value1[1]+=value2[1]; //size is sum of both
			value2[1]=value1[0]; //smaller component is pointing to the root of the larger component			
		}
	}
	
	public int find( T t ) {
		int root;
		int[] initial = uf.get(t);
		if( initial == null ) {
			root=uf.size();
			uf.put(t, new int[]{root,-1});
			index.add(t);
		} else {
			T help=t;
			while( initial[1] >= 0 ) {
				help = index.get(initial[1]);
				initial = uf.get(help);
			}
			root = initial[0];
			
			help=t;
			initial = uf.get(help);
			while( initial[1] >= 0 ) {
				help = index.get(initial[1]);
				initial[1]=root;
				initial = uf.get(help);
			}
		}
		return root;
	}
	
	public HashMap<Integer,ArrayList<T>> getComponents() {
		HashMap<Integer,ArrayList<T>> hash = new HashMap<Integer,ArrayList<T>>();
		for( T current: index ) {
			int comp = find(current);
			ArrayList<T> list = hash.get(comp);
			if( list == null ) {
				list = new ArrayList<T>();
				hash.put(comp, list);
			}
			list.add(current);
		}
		return hash;
	}
}
