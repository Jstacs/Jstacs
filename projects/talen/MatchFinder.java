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

import java.util.HashMap;
import java.util.Iterator;

import projects.tals.TALgetterDiffSM;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;

public abstract class MatchFinder {

	protected HashMap<HashEntry, LimitedSortedList<Match>> scoreHash;
	protected HashMap<HashEntry, LimitedSortedList<Match>> scoreHashRc;
	
	public MatchFinder(){
		this.scoreHash = new HashMap<HashEntry, LimitedSortedList<Match>>();
		this.scoreHashRc = new HashMap<HashEntry, LimitedSortedList<Match>>();
	}
	
	public synchronized void hash(HashEntry en, LimitedSortedList<Match> matches, boolean rc){
		if(rc){
			if(scoreHashRc.size() > 32){
				Iterator<HashEntry> it = scoreHashRc.keySet().iterator();
				int i=0;
				while(it.hasNext() && i < 20){
					it.next();
					it.remove();
				}
			}
			scoreHashRc.put( en, matches );
		}else{
			if(scoreHash.size() > 32){
				Iterator<HashEntry> it = scoreHash.keySet().iterator();
				int i=0;
				while(it.hasNext() && i < 20){
					it.next();
					it.remove();
				}
			}
			scoreHash.put( en, matches );
		}
	}
	
	public synchronized LimitedSortedList<Match> getHashed(HashEntry en, boolean rc){
		if(rc){
			return scoreHashRc.get( en );
		}else{
			return scoreHash.get( en );
		}
	}
	
	public abstract LimitedSortedList<Match> getScoresAbove( Sequence tal, double t, int cap, boolean capBest, boolean rc );

	public static MatchFinder getMatchFinder(DataSet ds, TALgetterDiffSM model, int maxLength){
		long n = 0;
		for(int i=0;i<ds.getNumberOfElements();i++){
			n += ds.getElementAt( i ).getLength();
		}
		if(n < 5E5){
			System.out.println("simple");
			return new SimpleMatchFinder( ds, model );
		}else if(n < 1E8){
			System.out.println("pst");
			return new PartialStringTree(ds,Math.min( 11, maxLength ),Math.min( 14, maxLength ),model);
		}else{
			System.out.println("infix");
			return new InfixMatchFinder(ds,Math.min( 8, maxLength ),model);
		}
	}
	
	public void reset() {
		scoreHash.clear();
		scoreHashRc.clear();
	}
	
	public static class HashEntry{
		
		private Sequence tal;
		private Double usedThresh;
		private Integer usedCap;
		private Boolean capBest;
		/**
		 * @param tal
		 * @param usedThresh
		 * @param usedCap
		 * @param capBest
		 */
		public HashEntry( Sequence tal, double usedThresh, int usedCap, boolean capBest ) {
			this.tal = tal;
			this.usedThresh = usedThresh;
			this.usedCap = usedCap;
			this.capBest = capBest;
		}
		
		public boolean equals(Object o){
			if(o instanceof HashEntry){
				HashEntry en = (HashEntry)o;
				return en.tal.equals( tal ) && en.usedThresh.equals( usedThresh ) && en.usedCap.equals( usedCap ) && en.capBest.equals( capBest );
			}
			return false;
		}
		
		public Sequence getTal() {
			return tal;
		}
		
		public double getUsedThresh() {
			return usedThresh;
		}
		
		public int getUsedCap() {
			return usedCap;
		}
		
		public boolean isCapBest() {
			return capBest;
		}

		@Override
		public int hashCode() {
			int hash = tal.hashCode();
			hash = hash*31 + usedThresh.hashCode();
			hash = hash*31 + usedCap.hashCode();
			hash = hash*31 + capBest.hashCode();
			return hash;
		}
		
		
		
	}
	
	public static class Match{
		
		private int seqIdx;
		private int seqPos;
		private boolean rc;
		private Sequence tal;
		
		
		/**
		 * @param seqIdx
		 * @param seqPos
		 * @param score
		 * @param rc
		 */
		public Match( int seqIdx, int seqPos,  boolean rc ) {
			this.seqIdx = seqIdx;
			this.seqPos = seqPos;
			this.rc = rc;
		}

		public int getSeqIdx() {
			return seqIdx;
		}
		
		public int getSeqPos() {
			return seqPos;
		}

		
		public boolean isRc() {
			return rc;
		}
		
		
		public String toString(){
			return "("+seqIdx+", "+seqPos+", "+rc+")";
		}

		public Sequence getTal() {
			return tal;
		}

		public void setTal(Sequence tal) {
			this.tal = tal;
		}
		
		
		
	}
	
}