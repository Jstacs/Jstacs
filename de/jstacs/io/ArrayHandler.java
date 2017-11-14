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

package de.jstacs.io;

import java.lang.reflect.Array;
import java.lang.reflect.Method;
import java.util.HashSet;
import java.util.LinkedList;

/**
 * This class handles arrays with elements of generic type and enables the user
 * to cast, clone, and create arrays easily.
 * 
 * @author Jens Keilwagen
 */
public final class ArrayHandler {

	/**
	 * This method returns the deepest class in the class hierarchy that is the
	 * class or a superclass of all instances in the array.
	 * 
	 * @param <T>
	 *            the type of the array elements
	 * @param o
	 *            the array
	 * 
	 * @return the superclass of all elements of the given array
	 */
	@SuppressWarnings( "unchecked" )
	public static <T> Class getSuperClassOf( T... o ) {
		Class current;
		LinkedList<Class> classHierarchy = new LinkedList<Class>();
		HashSet<Class> hash = new HashSet<Class>();
		for( int i = 0; i < o.length; i++ ) {
			if( o[i] != null ) {
				current = o[i].getClass();
				if( classHierarchy.size() == 0 ) {
					while( current != Object.class ) {
						classHierarchy.add( current );
						hash.add( current );
						current = current.getSuperclass();
					}
					classHierarchy.add( current );
					hash.add( current );
				} else {
					while( !hash.contains( current ) ) {
						current = current.getSuperclass();
					}
					while( classHierarchy.get( 0 ) != current ) {
						hash.remove( classHierarchy.remove( 0 ) );
					}
					if( classHierarchy.size() == 1 ) {
						break;
					}
				}
			}
		}
		if( classHierarchy.size() > 0 ) {
			return classHierarchy.get( 0 );
		} else {
			return o.getClass().getComponentType();
		}
	}
	
	/**
	 * This method creates a new array of the superclass of all elements of the
	 * given array and copies the elements. The order of the elements is not
	 * changed.
	 * 
	 * <br>
	 * <br>
	 * 
	 * Here is an example that demonstrates what can be done:<br>
	 * <br>
	 * <code>
	 * Object[] o = { new UniformTrainSM( alphabetContainer ), new UniformTrainSM( alphabetContainer ) };<br>
	 * AbstractTrainableStatisticalModel[] a = (AbstractTrainableStatisticalModel[]) ArrayHandler.cast( o );<br>
	 * </code> <br>
	 * This should work fine, while<br>
	 * 
	 * <code>AbstractTrainableStatisticalModel[] a = (AbstractTrainableStatisticalModel[]) o;</code><br>
	 * 
	 * will not.
	 * 
	 * @param o
	 *            the array
	 * 
	 * @return the casted array with the copied elements
	 * 
	 * @see ArrayHandler#getSuperClassOf(Object[])
	 * @see ArrayHandler#cast(Class, Object[])
	 */
	@SuppressWarnings( "unchecked" )
	public static Object[] cast( Object[] o ) {
		return cast( getSuperClassOf( o ), o );
	}

	/**
	 * This method returns an array of a user-specified class with all elements in the given array <code>o</code>.
	 * If the given array is already from the correct class it is directly returned, otherwise a new array is created
	 * and filled with the elements of <code>o</code>. In both cases, the order of the elements is not changed.
	 * 
	 * <br>
	 * <br>
	 * 
	 * Here is an example that demonstrates what can be done:<br>
	 * <br>
	 * <code>
	 * Object[] o = { new UniformTrainSM( alphabetContainer ), new UniformTrainSM( alphabetContainer ) };<br>
	 * AbstractTrainableStatisticalModel[] a = ArrayHandler.cast( AbstractTrainableStatisticalModel.class, o );<br>
	 * </code> <br>
	 * This should work fine, while<br>
	 * 
	 * <code>AbstractTrainableStatisticalModel[] a = (AbstractTrainableStatisticalModel[]) o;</code><br>
	 * 
	 * will not.
	 * 
	 * @param <S>
	 *            the type of the array elements
	 * @param c
	 *            the class of the array items
	 * @param o
	 *            the array
	 * 
	 * @return the casted array with the copied elements
	 * 
	 * @see ArrayHandler#getSuperClassOf(Object[])
	 */
	@SuppressWarnings( "unchecked" )
	public static <S> S[] cast( Class<? extends S> c, Object[] o ) {
		if( o.getClass().getComponentType().equals( c ) ) {
			return (S[]) o;
		} else {
			S[] res = (S[])Array.newInstance( c, o.length );
			for( int i = 0; i < o.length; i++ ) {
				res[i] = (S) o[i];
			}
			return res;
		}
	}

	/**
	 * This method returns a deep copy of a (multi-dimensional) array of {@link Cloneable}s or primitives. For arrays of {@link Cloneable}s, each element of the
	 * array is cloned using its <code>clone()</code>-method.
	 * 
	 * @param <T>
	 *            the type of the array elements
	 * @param t
	 *            the array
	 * 
	 * @return a deep copy of the given array
	 * 
	 * @throws CloneNotSupportedException
	 *             if an element could not be cloned
	 * 
	 * @see Cloneable
	 */
	@SuppressWarnings( "unchecked" )
	public static <T extends Cloneable> T[] clone( final T... t ) throws CloneNotSupportedException {
		return (T[]) deepClone( t );
	}
	
	@SuppressWarnings( "unchecked" )
	private static Object deepClone( final Object s ) throws CloneNotSupportedException{
		Object res = null;
		if( s != null ) {
			Class k = s.getClass();
			if( k.isArray() ){
				Class c = k.getComponentType();
				int l = Array.getLength( s );
				res = Array.newInstance( c, l );
				if( c.isPrimitive() ) {
					System.arraycopy( s, 0, res, 0, l );
				} else {
					for( int i = 0; i < l; i++ ) {
						Array.set( res, i, deepClone( Array.get( s, i ) ) );
					}
				}				
			} else {
				try {
					Method cloneMethod = k.getMethod( "clone" );
					res = cloneMethod.invoke( s );
				} catch ( Exception e ) {
					CloneNotSupportedException cnse = new CloneNotSupportedException( e.getMessage()+" for class "+k );
					cnse.setStackTrace( e.getStackTrace() );
					throw cnse;
				}
			}
		}
		return res;
	}
	
	/**
	 * This method returns an array of length <code>l</code> that has at each position a clone of <code>t</code>.
	 * 
	 * @param <T>
	 *            the type of the array elements that implements
	 *            {@link Cloneable}
	 * @param t
	 *            the original element
	 * @param l
	 * 			  the dimension of the new array
	 * 
	 * @return an array of length <code>l</code> that has at each position a clone of <code>t</code>
	 * 
	 * @throws CloneNotSupportedException
	 *             if <code>t</code> could not be cloned
	 * 
	 * @see Cloneable
	 */
	@SuppressWarnings( "unchecked" )
	public static <T extends Cloneable> T[] createArrayOf( T t, int l ) throws CloneNotSupportedException {
		Class c = t.getClass();
		T[] res = (T[])Array.newInstance( c, l );
		try {
			Method cloneMethod = c.getMethod( "clone" );
			for( int i = 0; i < l; i++ ) {
				res[i] = (T)cloneMethod.invoke( t );
			}
		} catch ( Exception e ) {
			CloneNotSupportedException cnse = new CloneNotSupportedException( e.getMessage() );
			cnse.setStackTrace( e.getStackTrace() );
			throw cnse;
		}
		return res;
	}
}
