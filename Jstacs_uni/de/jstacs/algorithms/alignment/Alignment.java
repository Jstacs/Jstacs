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

import de.jstacs.algorithms.alignment.cost.Costs;
import de.jstacs.algorithms.alignment.cost.Costs.Direction;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sequence;

/**
 * Class for gapped alignments using the Gotohs algorithm.
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
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
		 * The local version of an alignment aligns a subsequence of the first sequence with a subsequence of the second sequence.
		 */
		LOCAL;
	}
	
	private AlignmentType type;
	
	private int startS1, endS1, startS2, endS2;
	
	private Sequence s1, s2;

	private Costs costs;

	private D[][] d;

	private E[][] e;

	private F[][] f;
	
	/**
	 * The number of secondary diagonals at both sides of the diagonal.
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

		this.s1 = s1; this.startS1 = startS1; this.endS1 = endS1;
		this.s2 = s2; this.startS2 = startS2; this.endS2 = endS2;
		
		int l1 = endS1-startS1, l2 = endS2-startS2, start, end;
		
		//initialize & compute
		if( d == null || d.length < l1+1 ||d[0].length < l2+1 ) {
			d = new D[l1 + 1][l2 + 1];
			e = new E[l1 + 1][l2 + 1];
			f = new F[l1 + 1][l2 + 1];
			
			for( int i = 0; i <= l1; i++ ) {
				start = Math.max(0,i-offDiagonal);
				end = Math.min(l2,Math.max(l2,i+offDiagonal));
				for( int j = start; j <= end; j++ ) {
					e[i][j] = new E( i, j );
					f[i][j] = new F( i, j );
					d[i][j] = new D( i, j );
				}
			}
		} else {
			for( int i = 0; i <= l1; i++ ) {
				start = Math.max(0,i-offDiagonal);
				end = Math.min(l2,Math.max(l2,i+offDiagonal));
				for( int j = start; j <= end; j++ ) {
					e[i][j].compute();
					f[i][j].compute();
					d[i][j].compute();
				}
			}
		}

		//create results: aligned strings, ...
		Element curr = null;
		int endPos, startPos = 0;
		if( type == AlignmentType.LOCAL ) {
			endPos=-1;
			for( int i = 0; i <= l1; i++ ) {
				start = Math.max(0,i-offDiagonal);
				end = Math.min(l2,Math.max(l2,i+offDiagonal));
				for( int j = start; j <= end; j++ ) {
					if( curr == null || d[i][j].cost < curr.cost ) {
						curr = d[i][j];
						endPos=startS1+i;
					}
				}
			}
		} else {
			endPos=l1;
			
			if( e[l1][l2].cost < f[l1][l2].cost && e[l1][l2].cost < d[l1][l2].cost ) {
				curr = e[l1][l2];
			} else if( f[l1][l2].cost < d[l1][l2].cost ) {
				curr = f[l1][l2];
			} else {
				curr = d[l1][l2];
			}
		}

		double cost = curr.cost;

		StringBuffer b1 = new StringBuffer();
		StringBuffer b2 = new StringBuffer();

		AlphabetContainer cont = s1.getAlphabetContainer();

		while( curr.pre != null ) {
			if( curr.direction == Direction.DIAGONAL ) {
				b1.insert( 0, cont.getSymbol( startS1 + curr.i - 1, s1.discreteVal( startS1 + curr.i - 1 ) ) );
				b2.insert( 0, cont.getSymbol( startS2 + curr.j - 1, s2.discreteVal( startS2 + curr.j - 1 ) ) );
			} else if( curr.direction == Direction.LEFT ) {
				b1.insert( 0, '-' );
				b2.insert( 0, cont.getSymbol( startS2 + curr.j - 1, startS2 + s2.discreteVal( curr.j - 1 ) ) );
			} else if( curr.direction == Direction.TOP ) {
				b1.insert( 0, cont.getSymbol( startS1 + curr.j - 1, s1.discreteVal( startS1 + curr.i - 1 ) ) );
				b2.insert( 0, '-' );
			}
			if( curr.j == l2 && (curr.direction == Direction.DIAGONAL || curr.direction == Direction.LEFT) ) {
				endPos = startS1+curr.i;
			}
			if( (type==AlignmentType.LOCAL || curr.j == 1) && (curr.direction == Direction.DIAGONAL || curr.direction == Direction.TOP) ) {
				startPos = startS1+curr.i-1;
			}
			
			curr = curr.pre;
		}

		return new PairwiseStringAlignment( b1.toString(), b2.toString(), cost, startPos, endPos );
	}

	/**
	 * Class for a general element of the DP (<b>d</b>ynamic
	 * <b>p</b>rogramming)-matrices of the Needleman-Wunsch algorithm.
	 * 
	 * @author Jan Grau
	 */
	private abstract class Element {

		/**
		 * The row index of the element in the matrix.
		 */
		protected int i;

		/**
		 * The column index of the element in the matrix.
		 */
		protected int j;

		/**
		 * The costs of the alignment until this element, i.e.
		 * <code>s1(1..i),s2(1..j)</code>.
		 */
		double cost;

		/**
		 * The direction of the predecessor in the DP-matrix.
		 */
		protected Direction direction;

		/**
		 * The predecessor-element in the DP-matrix.
		 */
		protected Element pre;

		/**
		 * Creates a new element for row <code>i</code> and column
		 * <code>j</code>.
		 * 
		 * @param i
		 *            the row index for the new element
		 * @param j
		 *            the column index for the new element
		 */
		public Element( int i, int j ) {
			this.i = i;
			this.j = j;
			compute();
		}
		
		public abstract void compute();
	}

	//(mis-)match (== S in http://de.wikipedia.org/wiki/Gotoh-Algorithmus)
	private class D extends Element {

		private D( int i, int j ) {
			super( i, j );
		}
		
		public void compute() {
			if( i == 0 && j == 0 ) {
				cost = 0;
			} else if( i == 0 && j > 0 ) {
				pre = e[i][j];
				cost = pre.cost;
				direction = Direction.SELF;
			} else if( i > 0 && j == 0 ) {
				pre = f[i][j];
				cost = pre.cost;
				direction = Direction.SELF;
			} else {
				double diag = d[i - 1][j - 1].cost + costs.getCostFor( s1, s2, startS1+i, startS2+j, Direction.DIAGONAL );
				double left = e[i][j].cost;
				double top = f[i][j].cost;

				if( diag < left && diag < top ) {
					cost = diag;
					direction = Direction.DIAGONAL;
					pre = d[i - 1][j - 1];
				} else if( left < top ) {
					cost = left;
					direction = Direction.SELF;
					pre = e[i][j];
				} else {
					cost = top;
					direction = Direction.SELF;
					pre = f[i][j];
				}
				
				if( type == AlignmentType.LOCAL ) {
					if( 0 < cost ) {
						cost = 0;
						direction = Direction.SELF;
						pre = null;
					}
				}
			}
		}
	}

	//gap in first string (== H in http://de.wikipedia.org/wiki/Gotoh-Algorithmus)
	private class E extends Element {

		private E( int i, int j ) {
			super( i, j );
		}
		
		public void compute() {
			if( j == 0 || e[i][j-1] == null ) {
				cost = Double.POSITIVE_INFINITY;
			} else if( i == 0 && j > 0 ) {
				direction = Direction.LEFT;
				if( type==AlignmentType.LOCAL ) {
					cost = 0; 
					pre = null;
				} else {
					pre = e[i][j - 1];
					cost = costs.getGapCostsFor( j );					
				}
			} else {
				double elong = e[i][j - 1].cost + costs.getElongateCosts();
				double start = d[i][j - 1].cost + costs.getGapCostsFor( 1 );
				if( elong < start ) {
					cost = elong;
					direction = Direction.LEFT;
					pre = e[i][j - 1];
				} else {
					cost = start;
					direction = Direction.LEFT;
					pre = d[i][j - 1];
				}
			}
		}
	}

	//gap in second string (== V in http://de.wikipedia.org/wiki/Gotoh-Algorithmus)
	private class F extends Element {

		private F( int i, int j ) {
			super( i, j );
		}
		
		public void compute() {
			if( i == 0 || f[i-1][j] == null ) {
				cost = Double.POSITIVE_INFINITY;
			} else if( i > 0 && j == 0 ) {
				direction = Direction.TOP;
				cost = type!=AlignmentType.GLOBAL ? 0 : costs.getGapCostsFor( i );
				pre = type==AlignmentType.LOCAL ? null : f[i - 1][j];
			} else if ( type==AlignmentType.SEMI_GLOBAL && i > 0 && j == s2.getLength() ) {
				direction = Direction.TOP;
				pre = d[i - 1][j];
				cost = pre.cost;
			} else {
				double elong = f[i - 1][j].cost + costs.getElongateCosts();
				double start = d[i - 1][j].cost + costs.getGapCostsFor( 1 );
				if( elong < start ) {
					cost = elong;
					direction = Direction.TOP;
					pre = f[i - 1][j];
				} else {
					cost = start;
					direction = Direction.TOP;
					pre = d[i - 1][j];
				}
			}
		}
	}
}
