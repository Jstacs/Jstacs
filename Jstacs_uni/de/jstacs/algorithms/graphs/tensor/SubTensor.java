/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.algorithms.graphs.tensor;


/**
 * This Tensor can be used to extract or use only a part of a complete {@link Tensor}.
 * 
 * @author Jens Keilwagen
 */
public class SubTensor extends Tensor {

	private int offset;
	private Tensor t;
	
	/**
	 * This constructor creates a {@link SubTensor} using the {@link Tensor} <code>t</code> for the nodes <code> offset, offset+1, ..., offset+length-1</code>.
	 * 
	 * @param t the underlying {@link Tensor}
	 * @param offset the offset in the nodes, i.e., all nodes smaller than <code>offset</code> will not be used
	 * @param length the number of nodes which will be used starting at <code>offset</code>
	 */
	public SubTensor( Tensor t, int offset, int length ) {
		super( length, t.order );
		if( offset >= 0 ) {
			this.offset = offset; 
		} else {
			throw new IllegalArgumentException( "The offset has to be non-negative." );
		}
		if( t.getNumberOfNodes() < offset+length ){
			throw new IllegalArgumentException( "The length out of range." );
		}
		this.t = t;
	}
	
	@Override
	public int[] getMaximalEdgeFor( byte k, int child, int... parents ) {
		check( parents );
		check( child );
		modify( +1, parents );
		int[] res = t.getMaximalEdgeFor( k, child+offset, parents );
		modify( -1, parents );		
		return res;
	}

	@Override
	public double getRootValue( int child ) {
		check( child );
		return t.getRootValue( child+offset );
	}

	@Override
	public double getValue( byte k, int child, int... parents ) {
		check( parents );
		check( child );
		modify( +1, parents );
		double res = t.getValue( k, child+offset, parents );
		modify( -1, parents );
		return res;
	}

	@Override
	public void resetValue( byte k, int child, int... parents ) {
		check( parents );
		check( child );
		modify( +1, parents );
		t.resetValue( k, child+offset, parents );
		modify( -1, parents );		
	}

	@Override
	public void setRootValue( int child, double val ) {
		check( child );
		t.setRootValue( child+offset, val );
	}

	@Override
	public void setValue( byte k, double val, int child, int... parents ) {
		check( parents );
		check( child );
		modify( +1, parents );
		t.setValue( k, val, child+offset, parents );
		modify( -1, parents );
	}
	
	private void modify( int direction, int[] array ) {
		for( int i = 0; i < array.length; i++ )  {
			array[i] += direction*offset;
		}
	}
	
	private void check( int ... array ) throws IndexOutOfBoundsException {
		for( int i = 0; i < array.length; i++ )  {
			if( array[i] < 0 || array[i] >= getNumberOfNodes() ) {
				throw new IndexOutOfBoundsException( "check value: " + array[i] );
			}
		}
	}
}
