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

package de.jstacs.algorithms.alignment;

import de.jstacs.algorithms.alignment.cost.AffineCosts;
import de.jstacs.algorithms.alignment.cost.Costs;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.Sequence;

/**
 * Class for computing optimal alignments using Needleman-Wunsch algorithm of for affine gap costs Gotohs algorithm.
 * The class also provides the possibility of specifying the number of off-diagonals. In this case, the current
 * implementation helps to reduce the runtime but not the memory consumption.
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @see AffineCosts
 */
public class Alignment {

	/**
	 * @author Jens Keilwagen
	 */
	public enum AlignmentType {
		/**
		 * The global version of an alignment aligns two complete sequences.
		 */
		GLOBAL,
		/**
		 * The semi global version of an alignment allows gaps at the begin and at the end of the second sequence without any costs.
		 */
		SEMI_GLOBAL,
		/**
		 * The free shift alignment allows gaps at the begin and at the end of both sequence without any costs. Be aware that the match cost should be negative.
		 */
		FREE_SHIFT,
		/**
		 * The local version of an alignment aligns a subsequence of the first sequence with a subsequence of the second sequence.
		 */
		LOCAL;
	}
	
	AlignmentType type;
	
	private int startS1, startS2;
	
	private Sequence s1, s2;
	int l1, l2;

	private Costs costs;	
	private AffineCosts aCosts;
	
	private AlignmentAlgorithm algorithm;

	private double[][][] d;	
	/**
	 * The number of secondary diagonals at both sides of the main diagonal.
	 * If the alignment is performed for {@link Sequence}s of different length, 
	 * the number of secondary diagonals must be at least as large as the difference
	 * of the lengths.
	 * This allows to use a banded version of the alignment algorithm.
	 */
	private int offDiagonal;

	/**
	 * Creates a new {@link Alignment} instance that aligns the sequences
	 * <code>s1</code> and <code>s2</code> using the costs defined in
	 * <code>costs</code>.
	 * 
	 * @param type
	 *            the type of the alignment
	 * @param costs
	 *            the costs
	 */
	public Alignment( AlignmentType type, Costs costs ) {
		this(type,costs,Integer.MAX_VALUE);
	}
	
	/**
	 * Creates a new {@link Alignment} instance that aligns the sequences
	 * <code>s1</code> and <code>s2</code> using the costs defined in
	 * <code>costs</code> and a banded version of the alignment algorithm.
	 * 
	 * If the alignment is performed for {@link Sequence}s of different length, 
	 * the number of secondary diagonals must be at least as large as the difference
	 * of the lengths.
	 * 
	 * @param type
	 *            the type of the alignment
	 * @param costs
	 *            the costs
	 * @param offDiagonal
	 *            the number of secondary diagonals at both sides of the diagonal           
	 */
	public Alignment( AlignmentType type, Costs costs, int offDiagonal ) {
		this.type = type;
		this.costs = costs;
		if( costs instanceof AffineCosts ) {
			aCosts = (AffineCosts) costs;
			algorithm = new AffineAlignment();
		} else {
			aCosts = null;
			algorithm = new SimpleAlgorithm();
		}
		this.offDiagonal = offDiagonal;
	}

	/**
	 * Computes and returns the alignment of <code>s1</code> and <code>s2</code>
	 * ({@link #Alignment(AlignmentType, Costs)}).
	 * 
	 * @param s1 the first sequence 
	 * @param s2 the second sequence 
	 *  
	 * @return the alignment
	 */
	public PairwiseStringAlignment getAlignment( Sequence s1, Sequence s2 ) {
		return getAlignment( s1, 0, s1.getLength(), s2, 0, s2.getLength() );
	}
	
	/**
	 * Computes and returns the alignment of <code>s1</code> and <code>s2</code>
	 * ({@link #Alignment(AlignmentType, Costs)}).
	 * 
	 * @param s1 the first sequence 
	 * @param startS1 the start position in the first sequence
	 * @param endS1 the end position in the first sequence
	 * @param s2 the second sequence 
	 * @param startS2 the start position in the second sequence
	 * @param endS2 the end position in the second sequence 
	 *  
	 * @return the alignment
	 */
	public PairwiseStringAlignment getAlignment( Sequence s1, int startS1, int endS1, Sequence s2, int startS2, int endS2 ) {

		this.s1 = s1; this.startS1 = startS1;
		this.s2 = s2; this.startS2 = startS2;
		
		l1 = endS1-startS1;
		l2 = endS2-startS2;
		
		if( Math.abs(l1-l2) > offDiagonal ) throw new IllegalArgumentException("The number of secondary diagonals must be at least as large as the difference of the lengths, but is "+offDiagonal+" < "+Math.abs(l1-l2)+".");
		
		int start, end, h;
		
		//initialize
		if( d == null || d[0].length < l1+1 || d[0][0].length < l2+1 ) {
			d = new double[aCosts == null ? 1 : 3][l1+1][l2+1];
		}
		
		//compute
		for( int i = 0; i <= l1; i++ ) {
			start = Math.max(0,i-offDiagonal);
			h = i+offDiagonal;
			if( h > 0 ) {//due to possible overflow
				end = Math.min(l2,h);
			} else {
				end = l2;
			}
			algorithm.reset( i, 0, start );
			for( int j = start; j <= end; j++ ) {
				algorithm.compute(i,j);
			}
			if( end != l2 ) algorithm.reset( i, end+1, l2 );
		}

		//create results: aligned strings, ...
		int[] index = {-1,-1,-1};
		int endPos, startPos = 0;
		if( type == AlignmentType.LOCAL ) {
			endPos=-1;
			index[0] = 0;
			for( int i = 0; i <= l1; i++ ) {
				start = Math.max(0,i-offDiagonal);
				h = i+offDiagonal;
				if( h > 0 ) {//due to possible overflow
					end = Math.min(l2,h);
				} else {
					end = l2;
				}
				for( int j = start; j <= end; j++ ) {
					if( index[1] < 0 || d[0][i][j] < d[0][index[1]][index[2]] ) {
						index[1] = i;
						index[2] = j;
						endPos=startS1+i;
					}
				}
			}
		} else {
			endPos=l1;
			index[1] = l1;
			index[2] = l2;
			if( aCosts == null ) {
				index[0] = 0;
			} else {
				if( d[1][l1][l2] < d[2][l1][l2] && d[1][l1][l2] < d[0][l1][l2] ) {
					index[0] = 1;
				} else if( d[2][l1][l2] < d[0][l1][l2] ) {
					index[0] = 2;
				} else {
					index[0] = 0;
				}
			}
		}

		double cost = d[index[0]][index[1]][index[2]];

		StringBuffer b1 = new StringBuffer();
		StringBuffer b2 = new StringBuffer();

		AlphabetContainer cont = s1.getAlphabetContainer();

		int[] next = index.clone();
		while( true ) {
			algorithm.next(index, next);
			//System.out.println( Arrays.toString(index) + "\t" + Arrays.toString(next) );
			if( next[0] < 0 ) {
				break;
			}
			
			int row = next[1] - index[1];
			int column = next[2] - index[2];
			if( row == -1 && column == -1 ) {
				b1.insert( 0, cont.getSymbol( startS1 + index[1] - 1, s1.discreteVal( startS1 + index[1] - 1 ) ) );
				b2.insert( 0, cont.getSymbol( startS2 + index[2] - 1, s2.discreteVal( startS2 + index[2] - 1 ) ) );
			} else if( column == -1 ) {//&& row == 0;
				b1.insert( 0, '-' );
				b2.insert( 0, cont.getSymbol( startS2 + index[2] - 1, s2.discreteVal( startS2 + index[2] - 1 ) ) );
			} else if( row == -1 ) {
				b1.insert( 0, cont.getSymbol( startS1 + index[1] - 1, s1.discreteVal( startS1 + index[1] - 1 ) ) );
				b2.insert( 0, '-' );
			}
			if( index[2] == l2 && column == -1 ) {
				endPos = startS1+index[1];
			}
			if( (type==AlignmentType.LOCAL || index[2] == 1) && row == -1 ) {
				startPos = startS1+index[1]-1;
			}
			
			System.arraycopy( next, 0, index, 0, 3 );
		}
		return new PairwiseStringAlignment( b1.toString(), b2.toString(), cost, startPos, endPos );
	}
	
	private static interface AlignmentAlgorithm {
		public void compute( int i, int j );
		public void reset( int i, int startJ, int endJ );
		public void next( int[] current, int[] next );
	}
	
	private class SimpleAlgorithm implements AlignmentAlgorithm {
		public void compute( int i, int j ) {
			computeDirection( i, j );
		}
		
		public byte computeDirection( int i, int j ) {
			byte direction = -1;
			if( i == 0 && j == 0 ) {
				d[0][i][j] = 0;
			} else if( i == 0 && j > 0 ) {
				if( type != AlignmentType.LOCAL ) {
					direction = 1;
				}				
				d[0][i][j] = type != AlignmentType.GLOBAL ? 0 : d[0][i][j-1]+costs.getGapCosts();
			} else if( i > 0 && j == 0 ) {
				if( type != AlignmentType.LOCAL ) {
					direction = 2;
				}
				d[0][i][j] = type != AlignmentType.GLOBAL ? 0 :  d[0][i-1][j] + costs.getGapCosts();
			} else {
				double diag = d[0][i - 1][j - 1] + costs.getCostFor( s1, s2, startS1+i, startS2+j );
				double left = d[0][i][j - 1] + costs.getGapCosts();
				double top = d[0][i-1][j];
				if ( (i < l1 && j < l2) //inner part of the matrix 
						|| !(type==AlignmentType.SEMI_GLOBAL && j == l2)
						|| !(type==AlignmentType.FREE_SHIFT && (j==l2 || i == l1)) ) {
					top += costs.getGapCosts();
				}
				
				if( diag < left && diag < top ) {
					d[0][i][j] = diag;
					direction = 0;
				} else if( left < top ) {
					d[0][i][j] = left;
					direction = 1;
				} else {
					d[0][i][j] = top;
					direction = 2;
				}
			}
			
			if( type == AlignmentType.LOCAL && 0 < d[0][i][j] ) {
				d[0][i][j] = 0;
				direction = -1;
			}
			return direction;
		}

		@Override
		public void reset(int i, int startJ, int endJ) {
			if( startJ >= 0 && startJ < d[0][i].length ) {
				d[0][i][startJ] = Double.POSITIVE_INFINITY;
			}
			if( endJ > 0 && endJ <= d[0][i].length ) {
				d[0][i][endJ-1] = Double.POSITIVE_INFINITY;
			}
			//Arrays.fill( d[0][i], startJ, endJ, Double.POSITIVE_INFINITY );
		}

		@Override
		public void next( int[] index, int[] next ) {
			byte res = computeDirection( index[1], index[2] );
			if( res < 0 ) {
				next[0] = -1;
			} else {
				System.arraycopy(index, 0, next, 0, 3);
				switch( res ) {
					case 0: next[1]--; next[2]--;break;
					case 2: next[1]--; break;
					case 1: next[2]--; break;
				}
			}
		}		
	}
	
	private class AffineAlignment implements AlignmentAlgorithm {

		@Override
		public void compute(int i, int j) {
			computeGap1( i, j );
			computeGap2( i, j );
			computeMatchMisMatch( i, j );
		}
		
		public byte computeMatchMisMatch( int i, int j ) {
			byte direction = -99;
			if( i == 0 && j == 0 ) {
				d[0][i][j] = 0;
			} else if( i == 0 && j > 0 ) {
				d[0][i][j] = d[1][i][j];
				direction = 1;
			} else if( i > 0 && j == 0 ) {
				d[0][i][j] = d[2][i][j];
				direction = 2;
			} else {
				double diag = d[0][i - 1][j - 1] + costs.getCostFor( s1, s2, startS1+i, startS2+j );
				double left = d[1][i][j];
				double top = d[2][i][j];

				if( diag < left && diag < top ) {
					d[0][i][j] = diag;
					direction = 0;
				} else if( left < top ) {
					d[0][i][j] = left;
					direction = 1;
				} else {
					d[0][i][j] = top;
					direction = 2;
				}
			}
			if( type == AlignmentType.LOCAL && 0 < d[0][i][j] ) {
				d[0][i][j] = 0;
				direction = -1;
			}
			return direction;
		}

		public byte computeGap1( int i, int j ) {
			byte direction = -100;
			if( j == 0 ) {
				d[1][i][j] = Double.POSITIVE_INFINITY;
			} else if( i == 0 /*&& j > 0*/ ) {
				if( type==AlignmentType.LOCAL ) {
					direction = -1;
					d[1][i][j] = 0;
				} else {
					direction = 1;
					d[1][i][j] = type==AlignmentType.FREE_SHIFT ? 0 : aCosts.getGapCostsFor( j );					
				}
			} else if ( type==AlignmentType.FREE_SHIFT /*&& j > 0*/ && i == l1 ) {
				direction = 0;//TODO
				d[1][i][j] = d[0][i][j-1];
			} else {
				double elong = d[1][i][j - 1] + aCosts.getElongateCosts();
				double start = d[0][i][j - 1] + aCosts.getGapCostsFor( 1 );
				if( elong < start ) {
					d[1][i][j] = elong;
					direction = 1;
				} else {
					d[1][i][j] = start;
					direction = 0;
				}
			}
			return direction;
		}

		public byte computeGap2( int i, int j ) {
			byte direction = -101;
			if( i == 0 ) {
				d[2][i][j] = Double.POSITIVE_INFINITY;
			} else if( /*i > 0 &&*/ j == 0 ) {
				direction = (byte) (type==AlignmentType.LOCAL ? -1 : 2);
				d[2][i][j] = type!=AlignmentType.GLOBAL ? 0 : aCosts.getGapCostsFor( i );
			} else if ( (type==AlignmentType.SEMI_GLOBAL || type==AlignmentType.FREE_SHIFT) /*&& i > 0*/ && j == l2 ) {
				direction = 0;
				d[2][i][j] = d[0][i - 1][j];
			} else {
				double elong = d[2][i - 1][j] + aCosts.getElongateCosts();
				double start = d[0][i - 1][j] + aCosts.getGapCostsFor( 1 );
				if( elong < start ) {
					d[2][i][j] = elong;
					direction = 2;
				} else {
					d[2][i][j] = start;
					direction = 0;
				}
			}
			return direction;
		}

		@Override
		public void next(int[] index, int[] next) {
			byte b;
			switch( index[0] ) {
				case 0:
					b = computeMatchMisMatch( index[1], index[2] );
					if( b < 0 ) {
						next[0] = -1;
					} else {
						System.arraycopy(index, 0, next, 0, 3);
						if( b == 0 ) {
							next[1]--;
							next[2]--;
						} else {
							next[0] = b;
						}
					}
					break;
				case 1:
					b = computeGap1( index[1], index[2] );
					if( b < 0 ) {
						next[0] = -1;
					} else {
						System.arraycopy(index, 0, next, 0, 3);
						next[2]--;
						next[0]=b;
					}
					break;
				case 2:
					b = computeGap2( index[1], index[2] );
					if( b < 0 ) {
						next[0] = -1;
					} else {
						System.arraycopy(index, 0, next, 0, 3);
						next[1]--;
						next[0]=b;
					}
					break;
			}
		}

		@Override
		public void reset(int i, int startJ, int endJ) {
			for( int k = 0; k < d.length; k++ ) {
				if( startJ >= 0 && startJ < d[k][i].length ) {
					d[k][i][startJ] = Double.POSITIVE_INFINITY;
				}
				if( endJ > 0 && endJ <= d[k][i].length ) {
					d[k][i][endJ-1] = Double.POSITIVE_INFINITY;
				}
				//Arrays.fill( d[k][i], startJ, endJ, Double.POSITIVE_INFINITY );
			}
		}		
	}
}