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

package de.jstacs.io;

import java.lang.reflect.Array;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import de.jstacs.Singleton;
import de.jstacs.Singleton.SingletonHandler;
import de.jstacs.Storable;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;

/**
 * Class for parsing standard data types and arrays in and out of an XML
 * {@link java.io.File}. The methods with prefix <code>append</code> or
 * <code>add</code> are for encoding, while methods with prefix
 * <code>extract</code> are for decoding.
 * 
 * Supported types include:
 * <ul>
 * <li> {@link Boolean}</li>
 * <li> {@link Byte} </li>
 * <li> {@link Character} </li>
 * <li> {@link Class} </li>
 * <li> {@link Double} </li>
 * <li> {@link Enum} </li>
 * <li> {@link Float} </li>
 * <li> {@link Integer} </li>
 * <li> {@link Long} </li>
 * <li> {@link Short} </li>
 * <li> {@link Singleton} </li>
 * <li> {@link Storable} </li>
 * <li> {@link String} </li>
 * <li> arrays of primitives (<code>int</code>, <code>double</code>, <code>char</code>, ...)
 * <li> and any arrays containing only instances of these types. </li>
 * </ul>
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public final class XMLParser {
	
	private static final Class<StringBuffer> stringBufferClass = StringBuffer.class;
	private static final String ARRAY_TAG = "pos"; 
	private static final HashSet<Class<?>> simpleParsable;
	
	private static final String CLASS_NAME = "className";
	private static final String LENGTH = "length";
	private static final String VALUE = "val";
	private static final String NULL = "null";
	private static final String ENUM = "name";
	
	static {
		simpleParsable = new HashSet<Class<?>>();
		//primitive types
		simpleParsable.add( Boolean.TYPE );
		simpleParsable.add( Byte.TYPE );
		simpleParsable.add( Character.TYPE );
		simpleParsable.add( Double.TYPE );
		simpleParsable.add( Float.TYPE );
		simpleParsable.add( Integer.TYPE );
		simpleParsable.add( Long.TYPE );
		simpleParsable.add( Short.TYPE );
		//wrapper classes
		simpleParsable.add( Boolean.class );
		simpleParsable.add( Byte.class );
		simpleParsable.add( Character.class );
		simpleParsable.add( Double.class );
		simpleParsable.add( Float.class );
		simpleParsable.add( Integer.class );
		simpleParsable.add( Long.class );
		simpleParsable.add( Short.class );
		//String
		simpleParsable.add( String.class );
	}
	
	/**
	 * Frames the {@link StringBuffer} <code>source</code> with equal tags
	 * &quot;&lt; <code>tag</code>&gt;&quot; and &quot;&lt;<code>/tag</code>
	 * &gt;&quot;.
	 * 
	 * @param source
	 *            the source {@link StringBuffer} that should be encoded in XML
	 * @param tag
	 *            the tags by which the {@link StringBuffer} should be framed
	 * 
	 * @see #addTagsAndAttributes(StringBuffer, String, String)
	 */
	public static void addTags( StringBuffer source, String tag ) {
		addTagsAndAttributes( source, tag, null );
	}

	/**
	 * Frames the {@link StringBuffer} <code>source</code> with &quot;&lt;
	 * <code>tag attributes</code>&gt;&quot; and &quot;&lt;<code>/tag</code>
	 * &gt;&quot;.
	 * 
	 * @param source
	 *            the source {@link StringBuffer} that should be encoded in XML
	 * @param tag
	 *            the tags by which the {@link StringBuffer} should be framed
	 * @param attributes
	 *            <code>null</code> or some attributes, i.e. <code>value=&quot;100&quot; confidence=&quot;0&quot;</code>
	 *            
	 * @see #parseAttributes(Map)
	 */
	public static void addTagsAndAttributes( StringBuffer source, String tag, String attributes ) {
		addTagsAndAttributes(source, tag, attributes, -1);
	}
	
	public static void addTagsAndAttributes( StringBuffer source, String tag, String attributes, int indentation ) {
		int l = source.length(), i=indentation<0?0:indentation;
		String s = "<" + tag + (attributes==null?"":(" " + attributes)) +  ">" + (l>0?"\n":"");
		insertIndentation(source, indentation, 0);
		source.insert( i, s );
		if( l>0 ) {
			addIndentation(source, indentation);
		}
		source.append( "</" + tag + ">\n" );
	}

	/**
	 * Appends an {@link Object} with the tags to the {@link StringBuffer} <code>xml</code>.
	 * 
	 * @param xml
	 *            the source {@link StringBuffer} that should be encoded in XML
	 * @param s
	 *            the object that should be framed by the tags and appended to the source,
	 *            <code>s</code> might be <code>null</code>, {@link Storable}, {@link String}, primitive type,
	 *            wrapper class of a primitive type, or any array thereof
	 * @param tag
	 *            the tags by which the value should be framed
	 *            
	 * @see #parseAttributes(Map)
	 * @see XMLParser#appendObjectWithTagsAndAttributes(StringBuffer, Object, String, String, boolean)
	 */
	public static void appendObjectWithTags( StringBuffer xml, Object s, String tag ) {
		appendObjectWithTagsAndAttributes( xml, s, tag, null );
	}

	/**
	 * Appends an {@link Object} with the tags and attributes to the {@link StringBuffer} <code>xml</code>.
	 * 
	 * @param xml
	 *            the source {@link StringBuffer} that should be encoded in XML
	 * @param s
	 *            the object that should be framed by the tags and appended to the source,
	 *            <code>s</code> might be <code>null</code>, {@link Storable}, {@link String}, primitive type,
	 *            wrapper class of a primitive type, or any array thereof
	 * @param tag
	 *            the tags by which the value should be framed
	 * @param attributes
	 *            <code>null</code> or some attributes, e.g., <code>value=&quot;100&quot; confidence=&quot;0&quot;</code>
	 *            
	 * @see #parseAttributes(Map)           
	 * @see XMLParser#appendObjectWithTagsAndAttributes(StringBuffer, Object, String, String, boolean)
	 */
	public static void appendObjectWithTagsAndAttributes( StringBuffer xml, Object s, String tag, String attributes ) {
		appendObjectWithTagsAndAttributes( xml, s, tag, attributes, true );
	}

	/**
	 * Appends an {@link Object} with the tags and attributes to the {@link StringBuffer} <code>xml</code>.
	 * 
	 * @param xml
	 *            the source {@link StringBuffer} that should be encoded in XML
	 * @param s
	 *            the object that should be framed by the tags and appended to the source,
	 *            <code>s</code> might be <code>null</code>, {@link Storable}, {@link String}, primitive type,
	 *            wrapper class of a primitive type, or any array thereof
	 * @param tag
	 *            the tags by which the value should be framed
	 * @param attributes
	 *            <code>null</code> or some attributes, e.g., <code>value=&quot;100&quot; confidence=&quot;0&quot;</code>
	 * @param writeClassInfo
	 *            a boolean to enable to write class information (e.g. {@link Class#getSimpleName()}) of the {@link Object} <code>s</code> into the XML.
	 *            
	 * @see #parseAttributes(Map)           
	 * @see XMLParser#extractObjectAndAttributesForTags(StringBuffer, String, Map, Map, Class)
	 */
	@SuppressWarnings( "unchecked" )
	public static void appendObjectWithTagsAndAttributes( StringBuffer xml, Object s, String tag, String attributes, boolean writeClassInfo ) {
		appendObjectWithTagsAndAttributes(xml, s, tag, attributes, writeClassInfo, -1);
	}
	
	public static final char INDENTATION_CHAR='\t';
	
	public static void addIndentation( StringBuffer xml, int indentation ) {
		insertIndentation(xml,indentation,-1);
	}
	
	public static void insertIndentation( StringBuffer xml, int indentation, int pos ) {
		for( int i = 0; i < indentation; i++ ) {
			if( pos < 0 ) {
				xml.append(INDENTATION_CHAR);
			} else {
				xml.insert(pos,INDENTATION_CHAR);
			}
		}
	}
	
	public static int nextIndentation( int indentation ) {
		return indentation + (indentation<0 ? 0 : 1);
	}
	
	public static void appendObjectWithTagsAndAttributes( StringBuffer xml, Object s, String tag, String attributes, boolean writeClassInfo, int indentation ) {
		addIndentation(xml,indentation);
		int nextIndentation = nextIndentation(indentation);
		xml.append( "<" + tag + (attributes==null?"":(" " + attributes)) + ">" );
		if( s == null ) {
			xml.append( NULL );
		} else {
			Class<? extends Object> k = s.getClass();
			if( writeClassInfo ) {
				appendObjectWithTagsAndAttributes( xml, k.getName(), CLASS_NAME, null, false, nextIndentation );
			}
			if( k.isArray() ){
				int l = Array.getLength( s );
				appendObjectWithTagsAndAttributes( xml, l, LENGTH, null, false, nextIndentation );
				Class c = k.getComponentType();
				if( simpleParsable.contains( c ) ) {
					writeClassInfo = false;
				}
				for( int i = 0; i < l; i++ ) {
					appendObjectWithTagsAndAttributes( xml, Array.get( s, i ), ARRAY_TAG, VALUE + "=\"" + i + "\"", writeClassInfo, nextIndentation );
				}
			} else {
				if( k.isEnum() ) {
					appendObjectWithTagsAndAttributes( xml, ((Enum)s).name(), ENUM, null, false, nextIndentation );
				} else if( s instanceof Class ) {
					appendObjectWithTagsAndAttributes( xml, ((Class)s).getName(), ENUM/*TODO*/, null, false, nextIndentation );
				} else if( s instanceof Singleton ) {
				} else if( simpleParsable.contains( k ) ) {
					xml.append( s instanceof String ? escape((String) s) : s );
				} else if( Storable.class.isInstance( s ) ) {
					xml.append( ((Storable) s).toXML() );
				} else {
					throw new RuntimeException( "Instance of " + k + " can not be parsed to XML." );
				}
			}
		}
		if( xml.charAt(xml.length()-1)=='\n' ) addIndentation(xml,indentation);
		xml.append( "</" + tag + ">" );
	}
	
	private static String[] split( String attrib ) {//TODO return list?
		ArrayList<String> res = new ArrayList<String>();
		int idx, start = -1, next;
		while( (idx = attrib.indexOf( '=',start+1 )) >= 0 ) {
			//there is another attribute
			char c = attrib.charAt(idx+1);
			next = idx+1;
			do {
				next = attrib.indexOf(c,next+1);
				if(next < 0){
					next = attrib.length()-1;
				}
			}while( next+1 < attrib.length() && attrib.charAt(next+1) != ' ' );
			
			res.add( attrib.substring(start+1,next+1) );
			start = next+1;
		}
		return res.toArray(new String[0]);
	}
	
	/**
	 * Parses the XML-attributes given in <code>attrs</code> and returns a {@link Map} of attributed names and associated values (as {@link String}s).
	 * @param attrs the list of XML-attributes
	 * @return the {@link Map} of attribute names and values
	 * @throws NonParsableException if any of the attributes does not have a value
	 * 
	 * @see #parseAttributes(Map)
	 */
	private static Map<String, String> parseAttributes( String attrs ) throws NonParsableException {
		Map<String, String> map = new TreeMap<String, String>();
		String[] parts = split( attrs );//attrs.split( "(?<!=)\\s+(?!=)" );
		int vallength = 0;
		for( String part : parts ) {
			part = part.trim();
			if( part.length() > 0 ) {
				String[] keyVal = part.split( "=" );
				if( keyVal.length != 2 ) {
					throw new NonParsableException( "Malformed attributes: " + attrs );
				}
				keyVal[0] = keyVal[0].trim();
				keyVal[1] = keyVal[1].trim();
				vallength = keyVal[1].length();
				char c = keyVal[1].charAt( 0 );
				if( (c == '"' || c == '\'') && keyVal[1].charAt( vallength - 1 ) == c ) {
					keyVal[1] = keyVal[1].substring( 1, vallength - 1 );
				}
				map.put( keyVal[0], keyVal[1] );
			}
		}
		return map;
	}
	
	/**
	 * This method parses a map of attribute, i.e. key-value-pairs, to a String representation that is used
	 * for encoding XML.
	 * 
	 * @param map the {@link Map} of attribute names and values
	 * @return the {@link String} representation of the attributes
	 * 
	 * @see #parseAttributes(String)
	 */
	public static String parseAttributes( Map<String, String> map ) {
		Iterator<Entry<String,String>> it = map.entrySet().iterator();
		Entry<String,String> e;
		StringBuffer res = new StringBuffer();
		while( it.hasNext() ) {
			e = it.next();
			res.append( e.getKey() + "=\"" + e.getValue() + "\" " );
		}
		int l = res.length();
		if( l > 0 ) {
			l--;
		}
		return res.substring(0, l);
	}

	/**
	 * Tests whether the attributes given in <code>filterAttributes</code> are present in <code>myAttrs</code> and if the values of present attributes are equal.
	 * @param myAttrs the attributes to be tested
	 * @param filterAttributes the attributes and associated values that must be present in <code>myAttrs</code>
	 * @return true if all attributes and values in <code>filterAttributes</code> could be found in <code>myAttrs</code>, false otherwise
	 */
	private static boolean testFilter( Map<String, String> myAttrs, Map<String, String> filterAttributes ) {

		Set<String> keys = filterAttributes.keySet();
		String myVal = null, filterVal = null;

		for( String key : keys ) {
			myVal = myAttrs.get( key );
			filterVal = filterAttributes.get( key );
			if( myVal != null || filterVal != null ) {
				if( myVal == null ) {
					return false;
				}
				myVal = myVal.trim();
				filterVal = filterVal.trim();
				if( !myVal.equals( filterVal ) ) {
					return false;
				}
			}
		}
		return true;
	}
	
	/**
	 * Extracts the contents of <code>source</code> between <code>tag</code> start and end tags.
	 * @param source the XML-code containing start and end tag
	 * @param tag the tag (without angle brackets)
	 * @return the contents of start and end tags, without these tags, as a {@link StringBuffer}
	 * @throws NonParsableException if start or end tag could not be found
	 */
	public static StringBuffer extractForTag( StringBuffer source, String tag ) throws NonParsableException {
		return extractForTag( source, tag, null, null );
	}
	
	/**
	 * Extracts the contents of <code>source</code> between <code>tag</code> start and end tags.
	 * If <code>attributes</code> is not <code>null</code>, the attributes of the start tag are added to this {@link Map}.
	 * If <code>filterAttributes</code> is not <code>null</code>, the start tag is accepted only if its attributes and associated values contain those defined in <code>filterAttributed</code>.
	 * 
	 * @param source the XML-code containing start and end tag
	 * @param tag the tag (without angle brackets)
	 * @param attributes a {@link Map} for attributes and values, or <code>null</code> if no attributes should be parsed.
	 * @param filterAttributes a {@link Map} of attributes and associated values, which must be present in the attributes of the start tag, or <code>null</code> for no filtering
	 * 
	 * @return the contents of start and end tags, without these tags, as a {@link StringBuffer}, if the XML does not contain such a tagged entry, the method returns <code>null</code>
	 * 
	 * @throws NonParsableException if the XML is malformed
	 */
	public static StringBuffer extractForTag( StringBuffer source, String tag, Map<String, String> attributes,
			Map<String, String> filterAttributes ) throws NonParsableException {
		return extractForTag(source, tag, attributes, filterAttributes, null);
	}
	
	private static StringBuffer extractForTag( StringBuffer source, String tag, Map<String, String> attributes,
			Map<String, String> filterAttributes, int[] index ) throws NonParsableException {

		//find start tag
		int start = findOpeningTag( source, tag, attributes, filterAttributes, index==null?0:index[1] );
		if( start < 0 ) {
			return null;
		}
		int endOfStart = source.indexOf( ">", start + tag.length() ) + 1;
		
		//find end tag
		int pos = findClosingTag( source, tag, endOfStart );
		int closepos = source.indexOf( ">", pos + 1 );
		
		//prepare result
		StringBuffer result = new StringBuffer( source.substring( endOfStart, pos ) );
		if( index == null ) {
			source.delete( start, closepos+1 ); //same as: pos+taglength+3 );
		} else {
			index[0] = start;
			index[1] = closepos+1;
		}
		return result;
	}
	
	/**
	 * This method allows to check whether an XML contains a tagged entry.
	 * 
	 * @param source the XML-code containing start and end tag
	 * @param tag the tag (without angle brackets)
	 * @param attributes a {@link Map} for attributes and values, or <code>null</code> if no attributes should be parsed.
	 * @param filterAttributes a {@link Map} of attributes and associated values, which must be present in the attributes of the start tag, or <code>null</code> for no filtering
	 * 
	 * @return <code>true</code> if the XML contains such a tag, otherwise <code>false</code>
	 * 
	 * @throws NonParsableException if the XML is malformed
	 */
	public static boolean hasTag( StringBuffer source, String tag, Map<String, String> attributes,
			Map<String, String> filterAttributes ) throws NonParsableException {
		
		int start = findOpeningTag( source, tag, attributes, filterAttributes, 0 );
		if( start < 0 ) {
			return false;
		}
		int endOfStart = source.indexOf( ">", start + tag.length() ) + 1;
		
		findClosingTag( source, tag, endOfStart );
		return true;
	}
	
	/**
	 * This method skips a certain entry in the XML.
	 * 
	 * @param source the XML-code containing start and end tag
	 * @param tag the tag (without angle brackets)
	 * @param attributes a {@link Map} for attributes and values, or <code>null</code> if no attributes should be parsed.
	 * @param filterAttributes a {@link Map} of attributes and associated values, which must be present in the attributes of the start tag, or <code>null</code> for no filtering
	 * 
	 * @return -1 if the XML does not contain such a tagged entry, otherwise the index of the position after the corresponding end tag
	 * 
	 * @throws NonParsableException if the XML is malformed
	 */
	private static int skipTag( StringBuffer source, String tag, Map<String, String> attributes,
			Map<String, String> filterAttributes ) throws NonParsableException {
		
		int start = findOpeningTag( source, tag, attributes, filterAttributes, 0 );
		if( start < 0 ) {
			return -1;
		}
		int endOfStart = source.indexOf( ">", start + tag.length() ) + 1;
		
		int pos = findClosingTag( source, tag, endOfStart );
		int closepos = source.indexOf( ">", pos + 1 );
		
		return closepos+1;
	}
	
	/**
	 * This method find the index of the opening tag.
	 * 
	 * @param source the XML-code containing start and end tag
	 * @param tag the tag (without angle brackets)
	 * @param attributes a {@link Map} for attributes and values, or <code>null</code> if no attributes should be parsed.
	 * @param filterAttributes a {@link Map} of attributes and associated values, which must be present in the attributes of the start tag, or <code>null</code> for no filtering
	 * @param offset the offset for starting the search
	 * 
	 * @return -1 if no opening tag is found, otherwise the index of the the opening tag
	 * 
	 * @throws NonParsableException if it is not possible to find the closing &quot;&gt;&quot; for &quot;&lt;tag&quot; 
	 */
	private static int findOpeningTag( StringBuffer source, String tag, Map<String, String> attributes,
			Map<String, String> filterAttributes, int offset ) throws NonParsableException {
		if( source == null || source.length() == 0 ) {
			return -1;
		}
		
		//find position of first match tag (& attributes)
		int taglength = tag.length(), start, endOfStart = offset;
		Map<String, String> myMap = null;
		boolean found = false;
		
		String startTag = "<";
		String endTag = ">";
		String currentTag;
		int idx;
		do{
			start = source.indexOf( startTag, endOfStart );
			if( start < 0 ){
	 			return -1;
			} else if( source.charAt(start+1) == '/' ) {
				throw new NonParsableException( "Malformed XML: Found unexpected end tag: " + source.substring(start, source.indexOf(endTag, start+1)+1) + "." ); 
			} else {
				endOfStart = source.indexOf( endTag, start );
				
				idx = source.indexOf( " ", start );
				if( idx < 0 || idx > endOfStart ) {
					idx = endOfStart;
				}
			}
			endOfStart += 1;
			
			currentTag = source.substring( start+1, idx );
			if( currentTag.equals( tag ) ) {
				if( filterAttributes != null || attributes != null ) {
					myMap = parseAttributes( source.substring( start + taglength + 2, endOfStart - 1 ) );
				}
				if( filterAttributes != null ) {
					found = testFilter( myMap, filterAttributes );
				} else  {
					found = true;
				}
				if( found && attributes != null ) {
					attributes.putAll( myMap );
				}
				if( found ) {
					return start;
				}
			}
			
			endOfStart = findClosingTag( source, currentTag, endOfStart );
			endOfStart += 3 + currentTag.length();
		} while( true );
	}
	
	/**
	 * This method find the index of the closing tag.
	 * 
	 * @param source the XML-code containing start and end tag
	 * @param tag the tag (without angle brackets)
	 * @param offset the offset for starting the search
	 * 
	 * @return the index of the the closing tag
	 * 
	 * @throws NonParsableException if the XML is malformed (e.g. end tag not found, wrong end tag, ...) 
	 */
	private static int findClosingTag( StringBuffer source, String tag, int offset ) throws NonParsableException {	
		int counter = 1, idx, h;
		String endTag = "</" + tag + ">", startPrefix = "<" + tag;
		do {
			idx = source.indexOf( endTag, offset );
			if( idx < 0 ) {
				throw new NonParsableException( "Malformed XML: No end tag found for " + tag + "." );
			} else {
				counter--;
			}
			
			//check for start tags in the interval
			while( (h = source.indexOf( startPrefix, offset )) >= offset &&  h < idx ) {
				offset = h + startPrefix.length();
				char c = source.charAt(offset);
				if( c == ' ' || c == '>' ) {
					counter++;
				}
			}
			
			offset = idx + endTag.length();
		} while( counter != 0 );
		return idx;
	}

	/**
	 * Returns the parsed value between the tags.
	 * 
	 * @param xml the source {@link StringBuffer} that should be decoded from XML
	 * @param tag the tags between which the value shall be taken
	 * 
	 * @return the value between the tags
	 * 
	 * @throws NonParsableException  if the value could not be parsed
	 * 
	 * @see XMLParser#extractObjectAndAttributesForTags(StringBuffer, String, Map, Map, Class)
	 */
	public static Object extractObjectForTags( StringBuffer xml, String tag ) throws NonParsableException{
		return extractObjectAndAttributesForTags( xml, tag, null, null );		
	}
	
	/**
	 * Returns the parsed value between the tags.
	 * 
	 * @param <T> the type of the parsed object
	 * @param xml the source {@link StringBuffer} that should be decoded from XML
	 * @param tag the tags between which the value shall be taken
	 * @param k the class used to parse the value between the tag, if <code>null</code> the class will be inferred from the XML
	 * 
	 * @return the value between the tags
	 * 
	 * @throws NonParsableException  if the value could not be parsed
	 * 
	 * @see XMLParser#extractObjectAndAttributesForTags(StringBuffer, String, Map, Map, Class)
	 */
	public static<T> T extractObjectForTags( StringBuffer xml, String tag, Class<T> k ) throws NonParsableException{
		return extractObjectAndAttributesForTags( xml, tag, null, null, k );
	}
	
	/**
	 * Returns the parsed value between the tags.
	 * 
	 * @param xml the source {@link StringBuffer} that should be decoded from XML
	 * @param tag the tags between which the value shall be taken
	 * @param attributes a {@link Map} which can be used to obtain the attribute of the tag
	 * @param filterAttributes a {@link Map} which defines a filter for the tags
	 * 
	 * @return the value between the tags
	 * 
	 * @throws NonParsableException  if the value could not be parsed
	 * 
	 * @see XMLParser#extractObjectAndAttributesForTags(StringBuffer, String, Map, Map, Class)
	 */
	public static Object extractObjectAndAttributesForTags( StringBuffer xml, String tag, Map<String, String> attributes,
			Map<String, String> filterAttributes ) throws NonParsableException{
		return extractObjectAndAttributesForTags( xml, tag, attributes, filterAttributes, null );
	}
	
	/**
	 * Returns the parsed value between the tags.
	 * 
	 * @param <T> the type of the parsed object
	 * @param xml the source {@link StringBuffer} that should be decoded from XML
	 * @param tag the tags between which the value shall be taken
	 * @param attributes a {@link Map} which can be used to obtain the attribute of the tag
	 * @param filterAttributes a {@link Map} which defines a filter for the tags
	 * @param k the class used to parse the value between the tag, if <code>null</code> the class will be inferred from the XML
	 * 
	 * @return the value between the tags
	 * 
	 * @throws NonParsableException  if the value could not be parsed
	 * 
	 * @see XMLParser#appendObjectWithTagsAndAttributes(StringBuffer, Object, String, String, boolean)
	 */
	public static<T> T extractObjectAndAttributesForTags( StringBuffer xml, String tag, Map<String, String> attributes,
			Map<String, String> filterAttributes, Class<T> k ) throws NonParsableException{
		return extractObjectAndAttributesForTags( xml, tag, attributes, filterAttributes, k, null, null );
		
	}
	
	/**
	 * Returns the parsed value between the tags as an inner instance of the object <code>outerInstance</code>.
	 * 
	 * @param <T> the type of the parsed object
	 * @param <S> the type of the outer object
	 * @param xml the source {@link StringBuffer} that should be decoded from XML
	 * @param tag the tags between which the value shall be taken
	 * @param attributes a {@link Map} which can be used to obtain the attribute of the tag
	 * @param filterAttributes a {@link Map} which defines a filter for the tags
	 * @param k the class used to parse the value between the tag, if <code>null</code> the class will be inferred from the XML
	 * @param outerClass the class of the outer instance of the parsed object
	 * @param outerInstance the outer instance of the parsed object
	 * 
	 * @return the value between the tags
	 * 
	 * @throws NonParsableException  if the value could not be parsed
	 * 
	 * @see XMLParser#appendObjectWithTagsAndAttributes(StringBuffer, Object, String, String, boolean)
	 */
	public static<T, S> T extractObjectAndAttributesForTags( StringBuffer xml, String tag, Map<String, String> attributes,
			Map<String, String> filterAttributes, Class<T> k, Class<S> outerClass, S outerInstance ) throws NonParsableException{
		return extractObjectAndAttributesForTags( xml, tag, attributes, filterAttributes, k, outerClass, outerInstance, null );
	}
	
	@SuppressWarnings( "unchecked" )
	private static<T, S> T extractObjectAndAttributesForTags( StringBuffer xml, String tag, Map<String, String> attributes,
			Map<String, String> filterAttributes, Class<T> k, Class<S> outerClass, S outerInstance, int[] index ) throws NonParsableException{
		T erg;
		StringBuffer ex = extractForTag( xml, tag, attributes, filterAttributes, index );
		if( ex == null ) {
			throw new NonParsableException( "Could not find \"" + tag + "\"." );
		} else if( ex.toString().trim().equals( NULL ) ) {
			return null;
		}
		String className = null;
		int mod = k != null ? k.getModifiers() : 0;
		boolean infer = //k == null || k == Object.class || (!k.isArray() && !simpleParsable.contains( k ) && ( Modifier.isInterface( mod ) || Modifier.isAbstract( mod ) ) );
			k==null
			|| (
					!( k.isArray() || simpleParsable.contains(k) || Modifier.isFinal(mod) )
					&&
					( Modifier.isAbstract(mod) || Modifier.isInterface(mod) || k == Object.class )
				);
		if( infer ) {
			className = extractObjectAndAttributesForTags( ex, CLASS_NAME, null, null, String.class, null, null, null );
			try {
				k = (Class<T>) Class.forName( className );
			}catch( Exception e ) {
				throw getNonParsableException( "Class \""+className+"\" not found.", e );
			}
		}
		if( k.isArray() ){
			index = new int[2];
			
			Class c = k.getComponentType();
			int l = (int) extractObjectAndAttributesForTags( ex, LENGTH, null, null, Integer.TYPE, null, null, index );
			erg = (T) Array.newInstance( c, l );

			if( !simpleParsable.contains( c ) ) {
				c = infer?null:c;
			}
			
			Map<String, String> myFilterAttributes = new TreeMap<String, String>();			
			for( int i = 0; i < l; i++ ) {
				myFilterAttributes.clear();
				myFilterAttributes.put( VALUE, ""+i );
				Object o = extractObjectAndAttributesForTags( ex, ARRAY_TAG, null, myFilterAttributes, c, outerClass, outerInstance, index );
				Array.set( erg, i, o );
			}
		} else {		
			if( k.isEnum() ) {
				erg = (T) Enum.valueOf( (Class) k, extractObjectAndAttributesForTags( ex, ENUM, null, null, String.class, null, null, null ) );
			} else if( simpleParsable.contains( k ) ) {
				int offset=0;
				if( !infer ) {
					try {
						offset=skipTag(ex, CLASS_NAME, null, null);
					}catch(Exception e){
					};
				}
				erg = (T) cast( k, offset <= 0 ? ex.toString() : ex.substring( offset ) );
			} else if( Singleton.class.isAssignableFrom(k)) {
				try {
					erg = (T) SingletonHandler.getSingelton( (Class<Singleton>) k );
				} catch ( Exception e ) {
					throw getNonParsableException( "You must provide a static field SINGELTON.", e );
				}
			} else if( Storable.class.isAssignableFrom(k)) {
				//String old = ex.toString();
				try {
					if( outerClass != null ) {
						erg = k.getConstructor( outerClass, stringBufferClass ).newInstance( outerInstance, ex );
					} else {
						erg = k.getConstructor( stringBufferClass ).newInstance( ex );
					}
				} catch ( NoSuchMethodException e ) { 
					throw getNonParsableException( "You must provide a constructor " + k.getName() + "(StringBuffer).", e );
				} catch ( Exception e ) {
					throw getNonParsableException( "problem at " + k.getName() + ": " + e.getClass().getSimpleName() + ": " + e.getCause().toString()+"\n"+Arrays.toString( e.getCause().getStackTrace() ).replaceAll( "(,|\\[|\\])", "\n-- " )+"]]", e );
				}
			} else if( k == Class.class ) {
				String str = null;
				String temp = ex.toString();
				try {
					str = extractObjectAndAttributesForTags( ex, ENUM/*TODO*/, null, null, String.class, null, null, null );
					erg = (T) Class.forName( str );
				} catch ( ClassNotFoundException e ) {
					throw getNonParsableException( "Could not find the specified class: "+str+" ("+k+": "+temp+").", e );
				}
			} else {
				throw new NonParsableException( "Could not parse the object with tag \"" + tag + "\" and class \"" + k + "\"." );
			}
		}
		return erg;
	}
	
	
	/**
	 * Stores a set of {@link Sequence}s to XML, although {@link Sequence} does not implement {@link Storable}.
	 * Use for general data storage is discouraged and should be limited to those cases, where a {@link Sequence} needs
	 * to the stored within another object implementing {@link Storable}.
	 * @param xml the XML buffer
	 * @param tag the tag
	 * @param seqs the sequences
	 */
	public static void appendSequencesWithTags(StringBuffer xml, String tag, Sequence... seqs){
		String[] temp = null;
		if(seqs != null){
			temp = new String[seqs.length];
			for(int i=0;i<seqs.length;i++){
				if(seqs[i] != null){
					SequenceAnnotation[] anns = seqs[i].getAnnotation();
					AlphabetContainer alphabet = seqs[i].getAlphabetContainer();
					String seqstr = seqs[i].toString( alphabet.getDelim(), 0, seqs[i].getLength() );
					StringBuffer tb = new StringBuffer();
					XMLParser.appendObjectWithTags( tb, anns, "anns" );
					XMLParser.appendObjectWithTags( tb, alphabet, "alphabet" );
					XMLParser.appendObjectWithTags( tb, seqstr, "seqstr" );
					temp[i] = tb.toString();
				}
			}
		}
		XMLParser.appendObjectWithTags( xml, temp, tag );
	}
	
	/**
	 * Extracts a set of sequences from their XML representation.
	 * Use for general data storage is discouraged and should be limited to those cases, where a {@link Sequence} needs
	 * to the stored within another object implementing {@link Storable}.
	 * @param xml the XML buffer
	 * @param tag the tag
	 * @return the sequences
	 * @throws NonParsableException if the XML could not be parsed
	 * @throws WrongAlphabetException if the alphabet does not fit the sequence information
	 */
	public static Sequence[] extractSequencesWithTags(StringBuffer xml, String tag) throws NonParsableException, WrongAlphabetException{
		String[] temp = (String[])XMLParser.extractObjectForTags( xml, tag );
		if(temp != null){
			Sequence[] seqs = new Sequence[temp.length];
			for(int i=0;i<seqs.length;i++){
				if(temp[i] != null){
					StringBuffer tb = new StringBuffer( temp[i] );
					SequenceAnnotation[] anns = (SequenceAnnotation[])XMLParser.extractObjectForTags( tb, "anns" );
					AlphabetContainer alphabet = (AlphabetContainer)XMLParser.extractObjectForTags( tb, "alphabet" );
					String seqstr = (String)XMLParser.extractObjectForTags( tb, "seqstr" );
					seqs[i] = Sequence.create( alphabet, seqstr, alphabet.getDelim() );
				}
			}
			return seqs;
		}else{
			return null;
		}
	}
	
	//parsing String to instance of simpleParsable
	private static<T> Object cast( Class<T> c, String value ) throws NonParsableException {
		if( c == java.lang.String.class ) {
			return unescape(value);
		} else if( c == java.lang.Character.TYPE ) {
			return value.charAt( 0 );//XXX
		} else if( c == java.lang.Boolean.TYPE || c == Boolean.class ) {
			return Boolean.parseBoolean( value );
		} else if( c == java.lang.Byte.TYPE || c == Byte.class ) {
			return Byte.parseByte( value );
		} else if( c == java.lang.Short.TYPE || c == Short.class ) {
			return Short.parseShort( value );
		} else if( c == java.lang.Integer.TYPE || c == Integer.class ) {
			return Integer.parseInt( value );
		} else if( c == java.lang.Long.TYPE || c == Long.class ) {
			return Long.parseLong( value );
		} else if( c == java.lang.Float.TYPE || c == Float.class ) {
			return Float.parseFloat( value );
		} else if( c == java.lang.Double.TYPE || c == Double.class ) {
			return Double.parseDouble( value );
		} else {
			throw new NonParsableException( "Could not parse \"" + value + "\" to " + c + "." );
		}
	}
	
	/**
	 * The method returns a {@link NonParsableException} with a given error
	 * message for an {@link Exception}.
	 * 
	 * @param s
	 *            the error message of the {@link NonParsableException}
	 * @param e
	 *            the {@link Exception}
	 * 
	 * @return a {@link NonParsableException} with given error message
	 */
	private static final NonParsableException getNonParsableException( String s, Exception e ) {
		NonParsableException n = new NonParsableException( s );
		n.setStackTrace( e.getStackTrace() );
		return n;
	}
	
	/**
	 * This method parses the <code>original</code> {@link String} to <code>null</code> if <code>original</code> equals &quot;null&quot;.
	 * 
	 * @param original the original {@link String}
	 * 
	 * @return <code>null</code> or <code>original</code>
	 */
	public static String parseString( String original ) {
		return (original != null && original.equals( NULL )) ? null : original;
	}

	private static String[][] table = {
        	{ "&", "&amp;" },
			{ "\"", "&quot;" },
	        { "'", "&apos;" },
	        { ">", "&gt;" },
	        { "<", "&lt;" },
	};
	
	/**
	 * Escape special characters from a String allowing to wrap the String by e.g. XML tags.
	 * 
	 * @param original the original String
	 * 
	 * @return the escaped String
	 * 
	 * @see #unescape(String)
	 */
	public static String escape( String original ) {
		if( original != null ) {
			for( int i = 0; i < table.length; i++ ) {
				original = original.replaceAll( table[i][0], table[i][1] );
			}
		}
		return original;
	}
	
	/**
	 * Unescapes a String
	 * 
	 * @param escaped the escaped String. Opposite of the method {@link #escape(String)}.
	 * 
	 * @return the original String 
	 * 
	 * @see #escape(String)
	 */
	public static String unescape( String escaped ) {
		if( escaped != null ) {
			for( int i = table.length-1; i>=0; i-- ) {
				escaped = escaped.replaceAll( table[i][1], table[i][0] );
			}
		}
		return escaped;
	}
}
