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

import projects.tals.TALgetterDiffSM;
import de.jstacs.data.DataSet;
import de.jstacs.data.DiscreteSequenceEnumerator;
import de.jstacs.data.sequences.Sequence;


public class InfixMatchFinder extends MatchFinder implements Cloneable {

	private DataSet ds;
	private int infixLength;
	private int[] powers;
	private TALgetterDiffSM model;
	private HashMap<HashEntry, Object[]> preps;
	//public static long time = 0;
	
	public InfixMatchFinder(DataSet ds, int infixLength, TALgetterDiffSM model){
		this.ds = ds;
		this.infixLength = infixLength;
		powers = new int[infixLength];
		if((int) model.getAlphabetContainer().getAlphabetLengthAt( 0 ) != 4){
			throw new RuntimeException();
		}
		powers[0] = 1;
		for(int i=1;i<infixLength;i++){
			powers[i] = powers[i-1]*4;
		}
		this.model = model;
		try{
			this.model.fix();
		}catch(Exception e){
			throw new RuntimeException( e );
		}
		preps = new HashMap<HashEntry, Object[]>();
	}
	
	public InfixMatchFinder clone() throws CloneNotSupportedException {
		InfixMatchFinder clone = (InfixMatchFinder)super.clone();
		clone.model = model.clone();
		clone.scoreHash = (HashMap<HashEntry, LimitedSortedList<Match>>)scoreHash.clone();
		clone.scoreHashRc = (HashMap<HashEntry, LimitedSortedList<Match>>)scoreHashRc.clone();
		return clone;
	}
	
	public void setDataSet(DataSet ds){
		this.ds = ds;
		reset();
	}
	
	public synchronized Object[] getPreps(Sequence tal, double thresh) {
		
		//long time = System.currentTimeMillis();
		Object[] prep = preps.get( new HashEntry( tal, thresh, 0, false ) );
		if(prep == null){
		//	System.out.println("prep "+tal);
			prep = prepare( model, tal, infixLength, thresh );
			if(preps.size() > 6){
				preps.clear();
			}
			preps.put( new HashEntry( tal, thresh, 0, false ), prep );
		}
		//InfixMatchFinder.time += System.currentTimeMillis() - time;
		return prep;
	}
	
	@Override
	public LimitedSortedList<Match> getScoresAbove( Sequence tal, double thresh, int cap, boolean capBest, boolean rc ) {
		//Time time = new RealTime();
		HashEntry en = new HashEntry( tal, thresh, cap, capBest );
		LimitedSortedList<Match> list = null;
		list = getHashed( en, rc );
		if(list == null ){
			//System.out.println("computing matches");
		//	int length = infixLength;

		//	int pl = powers[length-1];

			Object[] prep = getPreps(tal, thresh);

			//System.out.println(t.getElapsedTime());

			boolean[] bools = (boolean[])prep[1];
			double[] scores = (double[])prep[0];
			boolean[] bools2 = (boolean[])prep[3];
			double[] scores2 = (double[])prep[2];
			double rest = (Double)prep[4];


			list = new LimitedSortedList<Match>( cap );
			//System.out.println("scanning");
		//	long n=0;
			//	t.reset();
		//	int order = Math.max(1,model.getOrder()+1);
		//	int length2 = Math.min( length, tal.getLength()+1-length+order );

			

		//	int off2 = length+length2-order;
		//	int ll2 = tal.getLength()-length-length2+order+1;
		//	boolean llb = ll2 > 0;

		//	double thresh2 = thresh-rest;


			for(int i=0;i<ds.getNumberOfElements();i++){
				//System.out.println(i);
				Sequence seq = ds.getElementAt( i );
				/*if(list.getLength() >= cap){
					if(!capBest){
						System.out.println("break");
						break;
					}
					if(list.getWorstScore() > thresh){
						thresh = list.getWorstScore();
			//			thresh2 = thresh-rest;
					}
				}*/
				if(rc){
					try{
						seq = seq.reverseComplement();
					}catch(Exception e){
						throw new RuntimeException();
					}
				}
				fillMatches( seq, i, tal, bools, bools2, scores, scores2, rest, thresh, rc, list, 0, seq.getLength()-tal.getLength()+1 );
			}
			
			hash( en, list, rc );
		}
		//System.out.println(time.getElapsedTime());
		return list;
	}

	public void fillMatches(Sequence seq, int i, Sequence tal, boolean[] bools, boolean[] bools2, double[] scores, double[] scores2, double rest, double thresh, boolean rc, LimitedSortedList<Match> list, int start, int end ){
		
		
		int length = infixLength;
		int order = model.getOrder();
		int length2 = Math.min( length, tal.getLength()+1-length+order );
		int ll2 = tal.getLength()-length-length2+order+1;
		boolean llb = ll2 > 0;
		int pl = powers[length-1];
		int pl2 = powers[length2-1];
		double thresh2 = thresh-rest;
		int off2 = length+length2-order;
		if(end-start>0){
			int idx = getIndex( seq, start, length );
			int idx3 = getIndex( seq, start+length-order-1, length2 );
			//System.out.println("new: ");
			//System.out.println(seq.getSubSequence( start, end-start+tal.getLength()-1 ));
			for(int j=start+1,k=start+length+length2-order-1,l=start+length;j<end;j++,k++,l++){
				idx3 = idx3/4 + seq.discreteVal( k )*pl2;//7s

				if(bools[idx] && bools2[idx3]){//20s

					double sc = scores[idx] + scores2[idx3];//10s

					if(sc >= thresh2){

						if(llb){
							sc += model.getPartialLogScoreFor( tal, seq, j-1, off2, ll2 );
						}
						if(sc >= thresh && list.checkInsert( sc )){
							if(rc){
								list.insert( sc, new Match(i,seq.getLength()-j-tal.getLength(),rc) );
							}else{
								list.insert( sc, new Match(i,j-1,rc) );
							}
							//dummy = sc;
						}
					}

				}


				idx = idx/4 + seq.discreteVal( l )*pl;//7s
			}
			//dummy = idx;
		}
	}
	
	private Object[] prepare(TALgetterDiffSM model, Sequence tal, int length, double thresh) {
		//System.out.println("preparing "+tal+" "+thresh+" "+this);
		double[] scs = new double[tal.getLength()+1]; 
		model.getBestPossibleScore( tal, scs );
		

		double rest = 0;
		
		for(int i=length;i<scs.length;i++){
				rest += scs[i];
		}
		int al = (int)model.getAlphabetContainer().getAlphabetLengthAt( 0 );
		
		
		boolean[] bools = new boolean[(int)Math.pow( al, length )];
		double[] scores = new double[bools.length];
		
		
		DiscreteSequenceEnumerator dse = new DiscreteSequenceEnumerator( model.getAlphabetContainer(), length, false );
		int i=0;
		while(dse.hasMoreElements()){
			Sequence seq = dse.nextElement();
			double sc = model.getPartialLogScoreFor( tal, seq, 0, 0, length );
			bools[i] = sc + rest >= thresh;
			scores[i] = sc;
			i++;
		}
		
		rest = 0;
		for(int j=0;j<length;j++){
			rest += scs[j];
		}

		int order = Math.max( model.getOrder(), 1);
		int length2 = Math.min( length, tal.getLength()+1-length+order );
		
		double rest2 = 0;
		for(int j=length+length2-order;j<scs.length;j++){
			rest += scs[j];
			rest2 += scs[j];
		}
		tal = tal.getSubSequence( length-order-1 );
		length = length2;
		
		boolean[] bools2 = new boolean[(int)Math.pow( al, length )];
		double[] scores2 = new double[bools2.length];
		dse = new DiscreteSequenceEnumerator( model.getAlphabetContainer(), length, false );
		
		i=0;
		while(dse.hasMoreElements()){
			Sequence seq = dse.nextElement();
			try{
				seq = Sequence.create( seq.getAlphabetContainer(), "T"+seq.toString() );
			}catch(Exception dnh){
				dnh.printStackTrace();
				throw new RuntimeException();
			}
			double sc = model.getPartialLogScoreFor( tal, seq, 0, order+1, length-order );
			bools2[i] = sc + rest >= thresh;
			scores2[i] = sc;
			i++;
		}
		return new Object[]{scores,bools,scores2,bools2,rest2};
		
	}
	
	private final int getIndex(Sequence seq, int off, int length){
		int idx = 0;
		for(int i=0;i<length;i++){
			idx += seq.discreteVal( off+i )*powers[i];
		}
		return idx;
	}
	
	
}
