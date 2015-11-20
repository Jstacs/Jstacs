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

import projects.tals.TALgetterDiffSM;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;


public class SimpleMatchFinder extends MatchFinder {

	private DataSet ds;
	private TALgetterDiffSM model;
	
	public SimpleMatchFinder(DataSet ds, TALgetterDiffSM model){
		this.ds = ds;
		this.model = model;
		try{
			this.model.fix();
		}catch(Exception e){
			throw new RuntimeException( e );
		}
	}

	public void setDataSet(DataSet ds){
		this.ds = ds;
		reset();
	}
	
	@Override
	public LimitedSortedList<Match> getScoresAbove( Sequence tal, double thresh, int cap, boolean capBest, boolean rc ) {
		LimitedSortedList<Match> list = null;
		HashEntry en = new HashEntry( tal, thresh, cap, capBest );
		list = null;//getHashed( en, rc );
		if(list == null ){
			list = new LimitedSortedList<Match>( cap );


			for(int i=0;i<ds.getNumberOfElements();i++){
				//System.out.println(i);
				Sequence seq = ds.getElementAt( i );
				if(list.getLength() >= cap){
					if(!capBest){
						System.out.println("break");
						break;
					}
					if(list.getWorstScore() > thresh){
						thresh = list.getWorstScore();
					}
				}
				if(rc){
					try{
						seq = seq.reverseComplement();
					}catch(Exception e){
						throw new RuntimeException();
					}
				}
				if(seq.getLength()>tal.getLength()){

					for(int j=0;j<seq.getLength()-tal.getLength();j++){
						double sc = model.getPartialLogScoreFor( tal, seq, j, 0, tal.getLength()+1 );
						if(sc >= thresh && list.checkInsert( sc )){
							if(rc){
								list.insert( sc, new Match(i,seq.getLength()-j-1-tal.getLength(),rc) );
							}else{
								list.insert( sc, new Match(i,j,rc) );
							}
						}
					}

				}
			}
			hash( en, list, rc );
		}
		return list;
	}

	
	
	
}
