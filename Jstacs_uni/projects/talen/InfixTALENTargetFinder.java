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

import javax.naming.OperationNotSupportedException;

import projects.talen.MatchFinder.Match;
import projects.tals.TALgetterDiffSM;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.utils.ComparableElement;


public class InfixTALENTargetFinder implements Cloneable{

	private InfixMatchFinder finder;
	private DataSet ds;
	private TALgetterDiffSM model;
	
	public InfixTALENTargetFinder(DataSet ds, TALgetterDiffSM model, int maxLength){
		this.finder = new InfixMatchFinder( ds, maxLength, model );
		this.ds = ds;
		this.model = model;
	}
	
	public InfixTALENTargetFinder(DataSet ds, TALgetterDiffSM model, InfixMatchFinder finder){
		this.finder = finder;
		this.ds = ds;
		this.model = model;
	}
	
	public InfixTALENTargetFinder clone() throws CloneNotSupportedException {
		InfixTALENTargetFinder clone = (InfixTALENTargetFinder)super.clone();
		clone.model = model.clone();
		clone.finder = finder.clone();
		return clone;
	}
	
	public void setDataSet(DataSet ds){
		this.ds = ds;
		( (InfixMatchFinder)this.finder ).setDataSet( ds );
		
	}
	
	public LimitedSortedList<TALENMatch> getTALENMatches(Sequence tal1, Sequence tal2, double totalThresh, double singleThresh1, double singleThresh2, int minDist, int maxDist, int limit, boolean relativeScores, boolean nTerm1, boolean nTerm2, boolean onlyhetero) throws OperationNotSupportedException{
		LimitedSortedList<TALENMatch> list = new LimitedSortedList<TALENMatch>( limit );
		if(nTerm1 && nTerm2){
			fillMatchesNTerm( tal1, tal2, totalThresh, singleThresh1, singleThresh2, minDist, maxDist, list, relativeScores );
			fillMatchesNTerm( tal2, tal1, totalThresh, singleThresh2, singleThresh1, minDist, maxDist, list, relativeScores );
			if(!onlyhetero){
				fillMatchesNTerm( tal1, tal1, totalThresh, singleThresh1, singleThresh1, minDist, maxDist, list, relativeScores );
				fillMatchesNTerm( tal2, tal2, totalThresh, singleThresh2, singleThresh2, minDist, maxDist, list, relativeScores );
			}
		}else if(nTerm1){
			fillMatchesCNTerm( tal2, tal1, totalThresh, singleThresh2, singleThresh1, minDist, maxDist, list, relativeScores, false );
			fillMatchesCNTerm( tal1, tal2, totalThresh, singleThresh1, singleThresh2, minDist, maxDist, list, relativeScores, true );
			if(!onlyhetero){
				fillMatchesNTerm( tal1, tal1, totalThresh, singleThresh1, singleThresh1, minDist, maxDist, list, relativeScores );
				fillMatches( tal2, tal2, totalThresh, singleThresh2, singleThresh2, minDist, maxDist, list, relativeScores );
			}
		}else if(nTerm2){
			fillMatchesCNTerm( tal1, tal2, totalThresh, singleThresh1, singleThresh2, minDist, maxDist, list, relativeScores, false );
			fillMatchesCNTerm( tal2, tal1, totalThresh, singleThresh2, singleThresh1, minDist, maxDist, list, relativeScores, true );
			if(!onlyhetero){
				fillMatchesNTerm( tal2, tal2, totalThresh, singleThresh2, singleThresh2, minDist, maxDist, list, relativeScores );
				fillMatches( tal1, tal1, totalThresh, singleThresh1, singleThresh1, minDist, maxDist, list, relativeScores );
			}
		}else{
			fillMatches( tal1, tal2, totalThresh, singleThresh1, singleThresh2, minDist, maxDist, list, relativeScores );
			fillMatches( tal2, tal1, totalThresh, singleThresh2, singleThresh1, minDist, maxDist, list, relativeScores );
			if(!onlyhetero){
				fillMatches( tal1, tal1, totalThresh, singleThresh1, singleThresh1, minDist, maxDist, list, relativeScores );
				fillMatches( tal2, tal2, totalThresh, singleThresh2, singleThresh2, minDist, maxDist, list, relativeScores );
			}
		}
		return list;
	}
	
	
	private void fillMatches(Sequence tal1, Sequence tal2, double totalThresh, double singleThresh1, double singleThresh2, int minDist, int maxDist, LimitedSortedList<TALENMatch> list, boolean relativeScores) throws OperationNotSupportedException{
		ComparableElement<Match, Double>[] list1Fwd = finder.getScoresAbove( tal1, singleThresh1*(relativeScores ? tal1.getLength()+1 : 1), -100000, true, false ).getSortedList();//TODO cap
		//ComparableElement<Match, Double>[] list2Fwd = finder.getScoresAbove( tal2, singleThresh2*(relativeScores ? tal2.getLength()+1 : 1), 1000000, true, true ).getSortedList();//TODO cap
		/*if(list1Fwd.length > 0){
			Match m = list1Fwd[0].getElement();

			Sequence sub3 = ds.getElementAt( m.getSeqIdx() ).getSubSequence( m.getSeqPos(), tal1.getLength()+1 );

			System.out.println(list1Fwd.length+" "+sub3+" "+model.getMatchString( tal1, sub3 )+" "+list1Fwd[0].getWeight());
		}*/


		Object[] prep = finder.getPreps(tal2,singleThresh2*(relativeScores ? tal2.getLength()+1 : 1) );

		LimitedSortedList<Match> list2 = new LimitedSortedList<MatchFinder.Match>( maxDist-minDist+tal2.getLength()+2 );
		
		//System.out.println(t.getElapsedTime());

		boolean[] bools = (boolean[])prep[1];
		double[] scores = (double[])prep[0];
		boolean[] bools2 = (boolean[])prep[3];
		double[] scores2 = (double[])prep[2];
		double rest = (Double)prep[4];
		
		for(int i=0;i<list1Fwd.length;i++){
			Match curr = list1Fwd[i].getElement();
			double currScore = list1Fwd[i].getWeight();
			
			int seq = curr.getSeqIdx();
			int pos = curr.getSeqPos();
			
			Sequence seq2 = ds.getElementAt( seq );
			if(pos+tal1.getLength()+1+minDist+tal2.getLength()+1 < seq2.getLength()){
				int length = maxDist-minDist+tal2.getLength()+1;
				if(length > seq2.getLength()-pos-tal1.getLength()-minDist-1){
					length = seq2.getLength()-pos-tal1.getLength()-minDist-1;
				}
				//Sequence sub = seq2.getSubSequence( pos+tal1.getLength()+1+minDist, length ).reverseComplement();
				//System.out.println("sub1: "+sub);
				//Sequence sub = seq2.reverseComplement().getSubSequence( seq2.getLength()-(pos+tal1.getLength()+1+minDist)-length, length );
				//System.out.println("sub2: "+sub);
				Sequence rc = seq2.reverseComplement();
				
				int tl = tal2.getLength();
				list2.clear();
				finder.fillMatches( rc, seq, tal2, bools, bools2, scores, scores2, rest, singleThresh2*(relativeScores ? tal2.getLength()+1 : 1), true, list2, seq2.getLength()-(pos+tal1.getLength()+1+minDist)-length, length-tl + (seq2.getLength()-(pos+tal1.getLength()+1+minDist)-length+1) );
				//System.out.println("l2: "+list2.getLength());
				for(int j=0;j<list2.getLength();j++){
					ComparableElement<Match, Double> m = list2.getElementAt( j );
					double score2 = m.getWeight();
					if(score2/(relativeScores ? tal2.getLength()+1 : 1) > singleThresh2 && currScore/(relativeScores ? tal1.getLength()+1 : 1) + score2/(relativeScores ? tal2.getLength()+1 : 1) > totalThresh){
						list.insert( currScore/(relativeScores ? tal1.getLength()+1 : 1) + score2/(relativeScores ? tal2.getLength()+1 : 1), new TALENMatch( tal1, tal2, curr,m.getElement(), (byte)0 ) );
					}
				}
				
				/*int n=0;
				System.out.println(rc.getSubSequence( seq2.getLength()-(pos+tal1.getLength()+1+minDist)-length, length-tl+tal2.getLength() ));
				for(int j=0,k=seq2.getLength()-(pos+tal1.getLength()+1+minDist)-length;j<length-tl;j++,k++){
					double score2 = model.getPartialLogScoreFor( tal2, rc, k, 0, tal2.getLength()+1 );
					if(score2/(relativeScores ? tal2.getLength()+1 : 1) > singleThresh2 && currScore/(relativeScores ? tal1.getLength()+1 : 1) + score2/(relativeScores ? tal2.getLength()+1 : 1) > totalThresh){
						
						/*System.out.println(seq2.getSubSequence( pos, tal1.getLength()+maxDist+tal2.getLength()+1 ));
						System.out.println(seq2.getSubSequence( pos,tal1.getLength()+1 ));
						System.out.println(model.getMatchString( tal1, seq2.getSubSequence( pos, tal1.getLength()+1 ) ));
						System.out.println(tal1);
						System.out.println(sub.getSubSequence( j, tal2.getLength()+1 ));
						System.out.println(model.getMatchString( tal2, sub.getSubSequence( j, tal2.getLength()+1 ) ));
						System.out.println(tal2);
						System.out.println(currScore+" "+score2+" > "+singleThresh1 + ", "+singleThresh2+", "+totalThresh);
						//System.out.println(seq);
						System.out.println("###############################");*//*
						n++;
						//list.insert( currScore/(relativeScores ? tal1.getLength()+1 : 1) + score2/(relativeScores ? tal2.getLength()+1 : 1), new TALENMatch( tal1, tal2, curr,new Match(seq,pos+tal1.getLength()+1+minDist+length-j-(tal2.getLength()+1),true) ) );//TODO position
					}
				}
				System.out.println("n: "+n);*/
			}
		}
	}
	
	private void fillMatchesNTerm(Sequence tal1, Sequence tal2, double totalThresh, double singleThresh1, double singleThresh2, int minDist, int maxDist, LimitedSortedList<TALENMatch> list, boolean relativeScores) throws OperationNotSupportedException{
		ComparableElement<Match, Double>[] list1Fwd = finder.getScoresAbove( tal1, singleThresh1*(relativeScores ? tal1.getLength()+1 : 1), -100000, true, false ).getSortedList();//TODO cap
		//ComparableElement<Match, Double>[] list2Fwd = finder.getScoresAbove( tal2, singleThresh2*(relativeScores ? tal2.getLength()+1 : 1), 1000000, true, true ).getSortedList();//TODO cap
		/*if(list1Fwd.length > 0){
			Match m = list1Fwd[0].getElement();

			Sequence sub3 = ds.getElementAt( m.getSeqIdx() ).getSubSequence( m.getSeqPos(), tal1.getLength()+1 );

			System.out.println(list1Fwd.length+" "+sub3+" "+model.getMatchString( tal1, sub3 )+" "+list1Fwd[0].getWeight());
		}*/


		Object[] prep = finder.getPreps(tal2,singleThresh2*(relativeScores ? tal2.getLength()+1 : 1) );

		LimitedSortedList<Match> list2 = new LimitedSortedList<MatchFinder.Match>( maxDist-minDist+tal2.getLength()+2 );
		
		//System.out.println(t.getElapsedTime());

		boolean[] bools = (boolean[])prep[1];
		double[] scores = (double[])prep[0];
		boolean[] bools2 = (boolean[])prep[3];
		double[] scores2 = (double[])prep[2];
		double rest = (Double)prep[4];
		
		for(int i=0;i<list1Fwd.length;i++){
			Match curr = list1Fwd[i].getElement();
			double currScore = list1Fwd[i].getWeight();
			
			int seq = curr.getSeqIdx();
			int pos = curr.getSeqPos();
			
			Sequence seq2 = ds.getElementAt( seq );
			if(pos - minDist - (tal2.getLength()+1) >= 0){
		//	if(pos+tal1.getLength()+1+minDist+tal2.getLength()+1 < seq2.getLength()){
//				int length = maxDist-minDist+tal2.getLength()+1;
				int start = pos - maxDist;
				int end = pos - minDist;
				//if(length > seq2.getLength()-pos-tal1.getLength()-minDist-1){
/*				if(length > seq2.getLength() - (tal2.getLength()+1+minDist )){
					length = seq2.getLength()-pos-tal1.getLength()-minDist-1;
				}
*/
				if(start < tal2.getLength()+1){
					start = tal2.getLength()+1;
				}
				int temp = start;
				start = seq2.getLength() - end;
				end = seq2.getLength() - temp+2;
				//Sequence sub = seq2.getSubSequence( pos+tal1.getLength()+1+minDist, length ).reverseComplement();
				//System.out.println("sub1: "+sub);
				//Sequence sub = seq2.reverseComplement().getSubSequence( seq2.getLength()-(pos+tal1.getLength()+1+minDist)-length, length );
				//System.out.println("sub2: "+sub);
				Sequence rc = seq2.reverseComplement();
				
				int tl = tal2.getLength();
				list2.clear();
				try{
					finder.fillMatches( rc, seq, tal2, bools, bools2, scores, scores2, rest, singleThresh2*(relativeScores ? tal2.getLength()+1 : 1), true, list2, start, end );
				}catch(Exception e){
					e.printStackTrace();
					System.out.println(start+" "+end+" "+rc.getLength()+"; "+pos+" "+minDist+" "+(pos - minDist - (tal2.getLength()+1)));
					throw new RuntimeException( e );
				}
				//System.out.println("l2: "+list2.getLength());
				for(int j=0;j<list2.getLength();j++){
					ComparableElement<Match, Double> m = list2.getElementAt( j );
					double score2 = m.getWeight();
					if(score2/(relativeScores ? tal2.getLength()+1 : 1) > singleThresh2 && currScore/(relativeScores ? tal1.getLength()+1 : 1) + score2/(relativeScores ? tal2.getLength()+1 : 1) > totalThresh){
						list.insert( currScore/(relativeScores ? tal1.getLength()+1 : 1) + score2/(relativeScores ? tal2.getLength()+1 : 1), new TALENMatch( tal1, tal2, curr,m.getElement(), (byte)3 ) );
					}
				}
				
				/*int n=0;
				System.out.println(rc.getSubSequence( seq2.getLength()-(pos+tal1.getLength()+1+minDist)-length, length-tl+tal2.getLength() ));
				for(int j=0,k=seq2.getLength()-(pos+tal1.getLength()+1+minDist)-length;j<length-tl;j++,k++){
					double score2 = model.getPartialLogScoreFor( tal2, rc, k, 0, tal2.getLength()+1 );
					if(score2/(relativeScores ? tal2.getLength()+1 : 1) > singleThresh2 && currScore/(relativeScores ? tal1.getLength()+1 : 1) + score2/(relativeScores ? tal2.getLength()+1 : 1) > totalThresh){
						
						/*System.out.println(seq2.getSubSequence( pos, tal1.getLength()+maxDist+tal2.getLength()+1 ));
						System.out.println(seq2.getSubSequence( pos,tal1.getLength()+1 ));
						System.out.println(model.getMatchString( tal1, seq2.getSubSequence( pos, tal1.getLength()+1 ) ));
						System.out.println(tal1);
						System.out.println(sub.getSubSequence( j, tal2.getLength()+1 ));
						System.out.println(model.getMatchString( tal2, sub.getSubSequence( j, tal2.getLength()+1 ) ));
						System.out.println(tal2);
						System.out.println(currScore+" "+score2+" > "+singleThresh1 + ", "+singleThresh2+", "+totalThresh);
						//System.out.println(seq);
						System.out.println("###############################");*//*
						n++;
						//list.insert( currScore/(relativeScores ? tal1.getLength()+1 : 1) + score2/(relativeScores ? tal2.getLength()+1 : 1), new TALENMatch( tal1, tal2, curr,new Match(seq,pos+tal1.getLength()+1+minDist+length-j-(tal2.getLength()+1),true) ) );//TODO position
					}
				}
				System.out.println("n: "+n);*/
			}
		}
	}
	
	private void fillMatchesCNTerm(Sequence tal1, Sequence tal2, double totalThresh, double singleThresh1, double singleThresh2, int minDist, int maxDist, LimitedSortedList<TALENMatch> list, boolean relativeScores, boolean rc) throws OperationNotSupportedException{
		
		ComparableElement<Match, Double>[] list1Fwd = finder.getScoresAbove( tal1, singleThresh1*(relativeScores ? tal1.getLength()+1 : 1), -100000, true, rc ).getSortedList();//TODO cap


		Object[] prep = finder.getPreps(tal2,singleThresh2*(relativeScores ? tal2.getLength()+1 : 1) );

		LimitedSortedList<Match> list2 = new LimitedSortedList<MatchFinder.Match>( maxDist-minDist+tal2.getLength()+2 );
		
		//System.out.println(t.getElapsedTime());

		boolean[] bools = (boolean[])prep[1];
		double[] scores = (double[])prep[0];
		boolean[] bools2 = (boolean[])prep[3];
		double[] scores2 = (double[])prep[2];
		double rest = (Double)prep[4];
		
		for(int i=0;i<list1Fwd.length;i++){
			Match curr = list1Fwd[i].getElement();
			double currScore = list1Fwd[i].getWeight();
			
			int seq = curr.getSeqIdx();
			int pos = curr.getSeqPos();
			
			Sequence seq2 = ds.getElementAt( seq );
			if(pos + minDist + tal1.getLength()+tal2.getLength()+2 < seq2.getLength()){
		//	if(pos+tal1.getLength()+1+minDist+tal2.getLength()+1 < seq2.getLength()){
//				int length = maxDist-minDist+tal2.getLength()+1;
				int start = pos+tal1.getLength()+1+minDist;
				int end = pos+tal1.getLength()+1+maxDist;

				if(end + tal2.getLength()+1 >= seq2.getLength()){
					end = seq2.getLength()-tal2.getLength()-1;
				}
				if(rc){
					seq2 = seq2.reverseComplement();
					int temp = start;
					start = seq2.getLength()-end-tal2.getLength()-1;
					end = seq2.getLength()-temp-tal2.getLength()-1;
				}
				
				int tl = tal2.getLength();
				list2.clear();
				try{
					finder.fillMatches( seq2, seq, tal2, bools, bools2, scores, scores2, rest, singleThresh2*(relativeScores ? tal2.getLength()+1 : 1), rc, list2, start, end );
				}catch(Exception e){
					e.printStackTrace();
					//System.out.println(start+" "+end+" "+rc.getLength()+"; "+pos+" "+minDist+" "+(pos - minDist - (tal2.getLength()+1)));
					throw new RuntimeException( e );
				}
				//System.out.println("l2: "+list2.getLength());
				for(int j=0;j<list2.getLength();j++){
					ComparableElement<Match, Double> m = list2.getElementAt( j );
					double score2 = m.getWeight();
					if(score2/(relativeScores ? tal2.getLength()+1 : 1) > singleThresh2 && currScore/(relativeScores ? tal1.getLength()+1 : 1) + score2/(relativeScores ? tal2.getLength()+1 : 1) > totalThresh){
						list.insert( currScore/(relativeScores ? tal1.getLength()+1 : 1) + score2/(relativeScores ? tal2.getLength()+1 : 1), new TALENMatch( tal1, tal2, curr,m.getElement(), rc ? (byte)1 : (byte)2 ) );
					}
				}
				
				/*int n=0;
				System.out.println(rc.getSubSequence( seq2.getLength()-(pos+tal1.getLength()+1+minDist)-length, length-tl+tal2.getLength() ));
				for(int j=0,k=seq2.getLength()-(pos+tal1.getLength()+1+minDist)-length;j<length-tl;j++,k++){
					double score2 = model.getPartialLogScoreFor( tal2, rc, k, 0, tal2.getLength()+1 );
					if(score2/(relativeScores ? tal2.getLength()+1 : 1) > singleThresh2 && currScore/(relativeScores ? tal1.getLength()+1 : 1) + score2/(relativeScores ? tal2.getLength()+1 : 1) > totalThresh){
						
						/*System.out.println(seq2.getSubSequence( pos, tal1.getLength()+maxDist+tal2.getLength()+1 ));
						System.out.println(seq2.getSubSequence( pos,tal1.getLength()+1 ));
						System.out.println(model.getMatchString( tal1, seq2.getSubSequence( pos, tal1.getLength()+1 ) ));
						System.out.println(tal1);
						System.out.println(sub.getSubSequence( j, tal2.getLength()+1 ));
						System.out.println(model.getMatchString( tal2, sub.getSubSequence( j, tal2.getLength()+1 ) ));
						System.out.println(tal2);
						System.out.println(currScore+" "+score2+" > "+singleThresh1 + ", "+singleThresh2+", "+totalThresh);
						//System.out.println(seq);
						System.out.println("###############################");*//*
						n++;
						//list.insert( currScore/(relativeScores ? tal1.getLength()+1 : 1) + score2/(relativeScores ? tal2.getLength()+1 : 1), new TALENMatch( tal1, tal2, curr,new Match(seq,pos+tal1.getLength()+1+minDist+length-j-(tal2.getLength()+1),true) ) );//TODO position
					}
				}
				System.out.println("n: "+n);*/
			}
		}
	}
	
	public void reset(){
		finder.reset();
	}
	
	public static class TALENMatch{
		
		private Sequence tal1;
		private Sequence tal2;
		private Match match1;
		private Match match2;
		private byte cat;//0: C/C, 1: C/N, 2: N/C, 3: NN
		
		public TALENMatch(Sequence tal1, Sequence tal2, Match match1, Match match2, byte cat){
			this.tal1 = tal1;
			this.tal2 = tal2;
			this.match1 = match1;
			this.match2 = match2;
			this.cat = cat;
		}

		public byte getCat(){
			return cat;
		}
		
		public Sequence getTal1() {
			return tal1;
		}

		
		public Sequence getTal2() {
			return tal2;
		}

		
		public Match getMatch1() {
			return match1;
		}

		
		public Match getMatch2() {
			return match2;
		}
		
		public String toString(){
			return match1+" "+match2;
		}
		
	}
	
}
