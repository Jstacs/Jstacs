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

package projects.talen;

import java.io.PrintWriter;
import java.util.HashMap;
import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import projects.tals.TALgetterDiffSM;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;


public class PartialStringTree extends MatchFinder {

	
	//public static PrintWriter pwr;
	private DataSet data;
	private int maxDepth;
	private int startDepth;
	private TALgetterDiffSM model;
	private Node root;
	private int n;
	
	public PartialStringTree(DataSet ds, int startDepth, int maxDepth, TALgetterDiffSM model) {
		this.data = ds;
		this.maxDepth = maxDepth;
		this.startDepth = startDepth;
		this.model = model;
		try{
			this.model.fix();
		}catch(Exception e){
			throw new RuntimeException( e );
		}
		
		construct();
	}
	
	
	private void construct() {

		root = new InnerNode( (int)data.getAlphabetContainer().getAlphabetLengthAt( 0 ));
		n=0;
		for(int i=0;i<data.getNumberOfElements();i++){
			//if(i%10 == 0){
				//System.out.println(i+"/"+data.getNumberOfElements());
				/*Set<Integer> keys = ArrayPool.map.keySet();
				Iterator<Integer> it = keys.iterator();
				while(it.hasNext()){
					int j=it.next();
					System.out.print(j+": "+ArrayPool.map.get( j ).size()+", ");
				}
				System.out.println();*/
			//}
			Sequence seq = data.getElementAt( i );
			for(int j=0;j<seq.getLength()-maxDepth+1;j++){
				root.add(seq,i,j,0,startDepth);
				n++;
			}
		}
		//System.out.println("constructed");
		//printN( new PrintWriter( "/Users/dev/Downloads/Ninitial.txt" ) );
		
		//Thread.sleep( 10000 );
		root.prune( 10 );
		System.gc();
		//printN( new PrintWriter( "/Users/dev/Downloads/Npruned.txt" ) );
		//System.out.println("pruned");
		//Thread.sleep( 10000 );
		root.expand( 20, maxDepth, 0 );
		System.gc();
		/*printN( new PrintWriter( "/Users/dev/Downloads/Nexpanded.txt" ) );
		System.out.println("expanded");
		Thread.sleep( 10000 );*/
		//System.out.println(root);
		
	}
	
	public void countChildren(int[] count){
		((InnerNode)root).countChildren( count );
	}
	
	/* (non-Javadoc)
	 * @see projects.talGA.MatchFinder#getScoresAbove(de.jstacs.data.sequences.Sequence, double, int, java.util.List, java.util.List, boolean)
	 */
	@Override
	public LimitedSortedList<Match> getScoresAbove(Sequence tal, double t, int cap, boolean capBest, boolean rc){
		//Time time = new RealTime();
		HashEntry en = new HashEntry( tal, t, cap, capBest );
		LimitedSortedList<Match> l = getHashed( en, rc );
		if(l == null ){
			l = new LimitedSortedList<Match>(cap);
			int[] currSeq = new int[tal.getLength()+1];
			double[] scs = new double[tal.getLength()+1];
			double sc = model.getBestPossibleScore( tal, scs );
			if(rc){
				for(int i=1;i<scs.length;i++){
					scs[i] += scs[i-1];
				}
				root.getScoresAboveRc(tal, t, currSeq, scs, 0, model, 0.0, l, cap, capBest);
			}else{
				for(int i=scs.length-2;i>=0;i--){
					scs[i] += scs[i+1];
				}
				root.getScoresAbove(tal, t, currSeq, scs, 0, model, 0.0, l, cap, capBest);
			}
			hash( en, l, rc );
		}
		//System.out.println(tal+" "+time.getElapsedTime());
		return l;
	}

	
	private static class ArrayPool{
		
		static HashMap<Integer, LinkedList<int[]>> map = new HashMap<Integer, LinkedList<int[]>>(); 
		
		static LinkedList<LinkedList<int[]>> pool = new LinkedList<LinkedList<int[]>>();
		
		public static int[] getArrayOfSize(int size, int[] old){
			LinkedList<int[]> l = map.get( size );
			int[] res = null;
			if(l != null && l.size() > 0){
				res = l.pop();
				if(l.size() == 0){
					LinkedList<int[]> el = map.remove( size );
					if(pool.size() < 20){
						pool.add( el );
					}
				}
			}else{
				res = new int[size];
			}
			System.arraycopy( old, 0, res, 0, old.length );
			LinkedList<int[]> l2 = map.get( old.length );
			if(l2 == null){
				if(pool.size() > 0){
					l2 = pool.pop();
				}else{
					l2 = new LinkedList<int[]>();
				}
				map.put( old.length, l2 );
			}
			if(l2.size() < 10){	
				l2.add( old );
			}
			
			return res;
		}
		
	}
	
	private abstract class Node{
		protected int n;
		
		protected abstract void add(Sequence seq, int seqIdx, int start, int depth, int length);

		public abstract void getScoresAbove( Sequence tal, double t, int[] currSeq, double[] bestSc, int depth, TALgetterDiffSM model, double currScore, LimitedSortedList<Match> list, int cap, boolean capBest );
		
		public abstract void getScoresAboveRc( Sequence tal, double t, int[] currSeq, double[] bestSc, int depth, TALgetterDiffSM model, double currScore, LimitedSortedList<Match> list, int cap, boolean capBest );

		public abstract void printN( PrintWriter wr );
		
		public abstract void prune(int minN);
		
		public abstract Node expand(int minN, int maxLength, int depth);
		
	}
	
	private class Leaf extends Node{

		private int[] seqIdxs;
		private int[] starts;
		
		private Leaf(){
			n=0;
		}
		
		public Leaf( InnerNode node ) {
			seqIdxs = new int[node.n];
			starts = new int[node.n];
			int n = node.fill(seqIdxs,starts,0);
			if(n != node.n){
				throw new RuntimeException();
			}
			this.n = n;
		}

		@Override
		protected void add( Sequence seq, int seqIdx, int start, int depth, int length ) {
			if(this.seqIdxs == null){
				this.seqIdxs = new int[1];
				this.starts = new int[1];
			}
			if(n+1==seqIdxs.length){
				seqIdxs = ArrayPool.getArrayOfSize( (int)Math.ceil(seqIdxs.length*1.1), seqIdxs );
				starts = ArrayPool.getArrayOfSize( (int)Math.ceil( starts.length*1.1), starts );
			}
			seqIdxs[n] = seqIdx;
			starts[n] = start;
			n++;
		}

		@Override
		public void getScoresAbove( Sequence tal, double t, int[] currSeq, double[] bestSc, int depth, TALgetterDiffSM model, double currScore,
				LimitedSortedList<Match> list, int cap, boolean capBest) {
			if(list.getLength()>=cap && !capBest){
				return;
			}
			for(int i=0;i<n;i++){
				//Sequence seq = data.getElementAt( seqIdxs[i] ).getSubSequence( starts[i] );
				Sequence full = data.getElementAt( seqIdxs[i] );
				if(full.getLength()-starts[i]<tal.getLength()+1){
					continue;
				}
				//System.out.println(seq);
				double temp = model.getPartialLogScoreFor( tal, full, starts[i], depth, tal.getLength()+1-depth );
				//System.out.println(temp);
				temp += currScore;
				//System.out.println(temp);
				if(temp > t){
					list.insert( temp, new Match(seqIdxs[i],starts[i],false) );
				}
			}
			//printN(pwr);
		}
		
		@Override
		public void getScoresAboveRc( Sequence tal, double t, int[] currSeq, double[] bestSc, int depth, TALgetterDiffSM model, double currScore,
				LimitedSortedList<Match> list, int cap, boolean capBest ) {
			if(list.getLength()>=cap && !capBest){
				return;
			}
			for(int i=0;i<n;i++){
				//Sequence seq = data.getElementAt( seqIdxs[i] ).getSubSequence( starts[i] );
				Sequence full = null;
				try{
					full = data.getElementAt( seqIdxs[i] ).reverseComplement();
				}catch(OperationNotSupportedException doesnothappen){
					throw new RuntimeException();
				}
				int start = full.getLength()-starts[i]-tal.getLength()-1;
				if(start < 0 || full.getLength()-start<tal.getLength()+1){
					continue;
				}
				int order = model.getOrder();
				//System.out.println(seq);
				double temp = model.getPartialLogScoreFor( tal, full, start, 0, tal.getLength()+2-depth+order-1 );//TODO
				//System.out.println(temp);
				temp += currScore;
				//System.out.println(temp);
				if(temp > t){
					list.insert( temp, new Match(seqIdxs[i],starts[i],true) );
				}
			}
			//printN(pwr);
		}
		
		public void printN(PrintWriter wr){
			wr.println(n);
		}

		@Override
		public void prune( int minN ) {
			
		}

		@Override
		public Node expand( int minN, int maxLength, int depth ) {
			if(this.n > minN){
				InnerNode temp = new InnerNode( (int)data.getAlphabetContainer().getAlphabetLengthAt( 0 ) );
				for(int i=0;i<this.n;i++){
					temp.add( data.getElementAt( seqIdxs[i] ), seqIdxs[i], starts[i], depth, maxLength );
				}
				temp.prune( minN );
				return temp;
			}else{
				return this;
			}
		}
		
	}
	
	private class InnerNode extends Node{
		
		private Node[] children;
		
		private InnerNode(int alphabetSize){
			this.children = new Node[alphabetSize];
			this.n = 0;
		}
		
		public int fill( int[] seqIdxs, int[] starts, int off ) {
			for(int i=0;i<children.length;i++){
				if(children[i] != null){
					if(children[i] instanceof InnerNode){
						off = ((InnerNode)children[i]).fill( seqIdxs, starts, off );
					}else{
						Leaf l = (Leaf)children[i];
						for(int j=0;j<l.n;j++){
							seqIdxs[off+j] = l.seqIdxs[j];
							starts[off+j] = l.starts[j];
						}
						off += l.n;
					}
				}
			}
			return off;
		}

		protected void add(Sequence seq, int seqIdx, int start, int depth, int length){
			
			int c = seq.discreteVal( start+depth );
			if(children[c] == null){
				if(depth+1==length){
					children[c] = new Leaf();
				}else{
					children[c] = new InnerNode(children.length);
				}
				
			}
			this.n++;
			children[c].add( seq, seqIdx, start, depth+1, length );
		}
		
		public void countChildren(int[] counts){
			int c = 0;
			for(int i=0;i<children.length;i++){
				if(children[i] != null){
					c++;
					if(children[i] instanceof InnerNode){
						((InnerNode)children[i]).countChildren(counts);
					}
				}
			}
			counts[c]++;
		}
		
		public void printN(PrintWriter wr){
			for(int i=0;i<children.length;i++){
				if(children[i] != null){
					children[i].printN(wr);
				}
			}
		}
		
		public String toString(){
			StringBuffer sb = new StringBuffer();
			for(int i=0;i<children.length;i++){
				if(children[i] != null){
					sb.append( i+": "+children[i].toString()+"\n" );
				}
			}
			/*if(sb.length() == 0){
				for(int i=0;i<seqIdxs.length();i++){
					sb.append( "("+seqIdxs.get( i )+","+starts.get( i )+") " );
				}
			}*/
			sb.append( "\n" );
			return sb.toString();
		}

		@Override
		public void getScoresAbove( Sequence tal, double t, int[] currSeq, double[] bestSc, int depth, TALgetterDiffSM model, double currScore, LimitedSortedList<Match> list, int cap, boolean capBest ) {
			if(list.getLength()>=cap && !capBest){
				return;
			}
			for(int i=0;i<children.length;i++){
				if(children[i] != null){
					currSeq[depth] = i;
					double tempScore = currScore + model.getPartialLogScoreFor( tal, currSeq, 0, depth, 1 );

					if(list.getLength() >= cap){
						if(t < list.getWorstScore()){
							t = list.getWorstScore();
						}
					}
					
					if( tempScore + bestSc[depth+1] >= t ){
						children[i].getScoresAbove( tal, t, currSeq, bestSc, depth+1, model, tempScore, list, cap, capBest );
					}

				}
			}
		}
		
		@Override
		public void getScoresAboveRc( Sequence tal, double t, int[] currSeq, double[] bestSc, int depth, TALgetterDiffSM model, double currScore, LimitedSortedList<Match> list, int cap, boolean capBest ) {
			if(list.getLength()>=cap && !capBest){
				return;
			}
			for(int i=0;i<children.length;i++){
				if(children[i] != null){
					
					if(list.getLength() >= cap){
						if(t < list.getWorstScore()){
							t = list.getWorstScore();
						}
					}
					
					currSeq[currSeq.length-depth-1] = children.length-i-1;
					int order = model.getOrder();
					double tempScore = 0;
					if(depth >= order){
						tempScore = currScore + model.getPartialLogScoreFor( tal, currSeq, 0, tal.getLength()-depth+order, 1 );//TODO
					}
					
					if( depth <= order || ( tempScore + bestSc[depth+1-order] >= t ) ){
						children[i].getScoresAboveRc( tal, t, currSeq, bestSc, depth+1, model, tempScore, list, cap, capBest );
					}
					
					
				}
			}
		}

		@Override
		public void prune( int minN ) {
			for(int i=0;i<children.length;i++){
				if(children[i] != null){
					if(children[i].n > minN){
						children[i].prune( minN );
					}else if(children[i] instanceof InnerNode){
						children[i] = new Leaf((InnerNode)children[i]);
					}
				}
			}
		}

		@Override
		public Node expand( int minN, int maxLength, int depth ) {
			for(int i=0;i<children.length;i++){
				if(children[i] != null && children[i].n > minN){
					children[i] = children[i].expand( minN, maxLength, depth+1 );
				}
			}
			return this;
		}
		
	}

	public int getNumberOfSequences() {
		return n;
	}


	public void printN( PrintWriter printWriter ) {
		root.printN( printWriter );
		printWriter.close();
	}
	
}
