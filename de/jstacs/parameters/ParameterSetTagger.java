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

package de.jstacs.parameters;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map.Entry;

import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.utils.ComparableElement;

/**
 * This class implements a tagger for {@link Parameter} of {@link ParameterSet}. This enable a parameter usage in command line
 * tool. Argument can be parsed from the specific syntax <code>tag&lt;delimiter&gt;value</code>.
 * 
 * <br><br>
 * 
 * <b>This class accesses the {@link ParameterSet} via {@link ParameterSet#getParameterAt(int)} only once during creation of
 * a new instance (in the constructor). If the number of parameters changes after the creation of a specific instance, this
 * instance will not detect this.</b>
 * 
 * @author Jens Keilwagen
 */
public class ParameterSetTagger {

	private Hashtable<String, ComparableElement<Parameter,Integer>> tagParameterHashtable;
	
	/**
	 * The constructor creates an new instance by collecting and tagging all parameters of the {@link ParameterSet}s.
	 *   
	 * @param tags the unambiguous tags for all parameters
	 * @param sets the sets of {@link Parameter}s
	 * 
	 * @see ParameterSet#getNumberOfParameters()
	 * @see ParameterSet#getParameterAt(int)
	 */
	public ParameterSetTagger( String[] tags, ParameterSet... sets ) {
		tagParameterHashtable = new Hashtable<String, ComparableElement<Parameter,Integer>>();
		int t = 0, n, s = 0;
		for( ; s < sets.length; s++ ) {
			for( n = 0; n < sets[s].getNumberOfParameters(); n++, t++ ) {
				if( tagParameterHashtable.containsKey( tags[t] ) ) {
					throw new IllegalArgumentException( "The tags have to be unambiguous. See tag: " + tags[t] );
				} else {
					tagParameterHashtable.put( tags[t], new ComparableElement<Parameter, Integer>( sets[s].getParameterAt( n ), t ) );
				}
			}
		}
		if( tags.length != t ) {
			throw new IllegalArgumentException( "Not all tags have been used." );
		}
	}
	
	/**
	 * 
	 * @param delimiter the delimiter between each tag and value
	 * @param args the arguments, each argument has the form <code>tag&lt;delimiter&gt;value</code>
	 * 
	 * @throws IllegalArgumentException if any argument could not be assigned to a parameter
	 * @throws IllegalValueException if any argument could not be parsed
	 */
	public void fillParameters( String delimiter, String... args ) throws IllegalArgumentException, IllegalValueException {
		for( int idx, i = 0; i < args.length; i++ ) {
			idx = args[i].indexOf( delimiter );
			if( idx < 0 ) {
				throw new IllegalArgumentException( "Unable to assign argument " + i + " to a parameter: " + args[i] );
			} else {
				setValueFromTag( args[i].substring( 0, idx ), args[i].substring( idx+delimiter.length() ) );
			}
		}
	}
	
	/**
	 * This method returns the {@link Parameter} specified by the <code>tag</code>
	 * 
	 * @param tag the tag of the {@link Parameter}
	 * 
	 * @return the {@link Parameter} specified by the <code>tag</code>
	 * 
	 * @see ParameterSet#getParameterAt(int)
	 */
	public Parameter getParameterFromTag( String tag ) {
		return tagParameterHashtable.get( tag ).getElement();
	}
	
	/**
	 * @param tag the tag of the {@link Parameter}
	 * 
	 * @return <code>true</code> if the the parameter is set
	 * 
	 * @see #getParameterFromTag(String)
	 * @see Parameter#isSet()
	 */
	public boolean isSet( String tag ) {
		return getParameterFromTag( tag ).isSet();
	}
	
	/**
	 * This method returns the value of the {@link Parameter} specified by the <code>tag</code>.
	 * 
	 * @param tag the tag of the {@link Parameter}
	 * 
	 * @return the value of the {@link Parameter} specified by the <code>tag</code>
	 * 
	 * @see #getParameterFromTag(String)
	 * @see Parameter#getValue()
	 */
	public Object getValueFromTag( String tag ) {
		return getParameterFromTag( tag ).getValue();
	}
	
	/**
	 * This method returns the casted value of the {@link Parameter} specified by the <code>tag</code>.
	 * 
	 * @param tag the tag of the {@link Parameter}
	 * @param c the class that is used for casting 
	 * 
	 * @return the value of the {@link Parameter} specified by the <code>tag</code>
	 * 
	 * @param <T> the type of the class
	 * 
	 * @see #getValueFromTag(String)
	 */
	public <T> T getValueFromTag( String tag, Class<T> c ) {
		return (T) getValueFromTag( tag );
	}
	
	/**
	 * This method allows to easily set the value of a parameter defined by the tag.
	 * 
	 * @param tag the tag of the {@link Parameter}
	 * @param value the value to be set
	 * 
	 * @throws IllegalValueException if the value could not be set
	 * 
	 * @see #getParameterFromTag(String)
	 * @see Parameter#setValue(Object)
	 */
	public void setValueFromTag( String tag, Object value ) throws IllegalValueException {
		getParameterFromTag( tag ).setValue( value );
	}
	
	/**
	 * This method allows to check whether all tagged parameters that require a value also have some value.
	 * 
	 * @return <code>true</code> if each {@link Parameter} that requires a value also has a default value or has been set 
	 * 
	 * @see Parameter#hasDefaultOrIsSet()
	 */
	public boolean hasDefaultOrIsSet() {
		Iterator<Entry<String,ComparableElement<Parameter, Integer>>> it = tagParameterHashtable.entrySet().iterator();
		Parameter current = null;
		while( it.hasNext() ) {
			current = it.next().getValue().getElement();
			if( !current.isRequired() || current.hasDefaultOrIsSet() ) {
				current = null;
			} else {
				break;
			}
		}
		return current == null;
	}
	
	/**
	 * This method allows to get a String representation where the tagged parameters are sorted in some specific way.
	 * 
	 * @param ec the comparator for sorting
	 * 
	 * @return a {@link String} representation of this instance
	 * 
	 * @see KeyEntryComparator
	 * @see RankEntryComparator
	 */
	public String toString( Comparator<Entry<String,ComparableElement<Parameter, Integer>>>  ec ) {
		Iterator<Entry<String,ComparableElement<Parameter, Integer>>> it = tagParameterHashtable.entrySet().iterator();
		Entry<String,ComparableElement<Parameter, Integer>>[] all = new Entry[tagParameterHashtable.size()];
		int i=0;
		while( it.hasNext() ) {
			all[i] = it.next();
			i++;
		}
		Arrays.sort( all, ec );
		StringBuffer sb = new StringBuffer();
		for( i = 0; i < all.length; i++ ) {
			sb.append( all[i].getKey() + "\t... " + all[i].getValue().getElement().toString() + "\n" );
		}
		return sb.toString();
	}
	
	public String toString() {
		return toString( new RankEntryComparator<String,Parameter>() );
	}

	/**
	 * This class implements a comparator on {@link java.util.Map.Entry} that sorts by the key of the {@link java.util.Map.Entry}.
	 * This class allows, for instance, to sort the {@link Parameter} in a {@link ParameterSetTagger} by their tags.
	 * 
	 * @author Jens Keilwagen
	 *
	 * @param <K> the type of the key
	 * @param <V> the type of the value
	 */
	public static class KeyEntryComparator<K extends Comparable<K>,V> implements Comparator<Entry<K,V>> {
		
		/**
		 * Creates a {@link Comparator} that only compares the keys of {@link Entry}s.
		 */
		public KeyEntryComparator(){}; 
		
		public int compare( Entry<K, V> o1, Entry<K, V> o2 ) {
			return o1.getKey().compareTo( o2.getKey() );
		}		
	}
	
	/**
	 * This class implements a comparator on {@link java.util.Map.Entry} where value is a {@link ComparableElement} with weight {@link Integer}.
	 * This class allows, for instance, to sort the {@link Parameter} in a {@link ParameterSetTagger} as given in the constructor. 
	 * 
	 * @author Jens Keilwagen
	 *
	 * @param <K> the type of the key
	 * @param <V> the type of the value
	 */
	public class RankEntryComparator<K,V> implements Comparator<Entry<K,ComparableElement<V,Integer>>> {
		
		/**
		 * Creates a {@link Comparator} that only compares the ranks of {@link Entry}s where the rank is determined by the {@link ParameterSet#parameters}.
		 */
		public RankEntryComparator(){}; 
		
		public int compare( Entry<K,ComparableElement<V,Integer>> o1, Entry<K,ComparableElement<V,Integer>>o2 ) {
			return o1.getValue().getWeight() - o2.getValue().getWeight();
		}
	}	
}
