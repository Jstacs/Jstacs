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
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.parameters;

import java.util.Arrays;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor;
import de.jstacs.utils.IntList;

/**
 * Class for a {@link Parameter} that provides a collection of possible values.
 * The user can select one or more values out of the collection as values of
 * this {@link Parameter}.
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @see AbstractSelectionParameter
 */
public class MultiSelectionParameter extends AbstractSelectionParameter {

	private boolean[] selected;
	private boolean[] defaultSelected;
	private int current;
	
	/**
	 * Constructor for a {@link MultiSelectionParameter}. The first
	 * option in the selection is selected by default.
	 * 
	 * @param datatype
	 *            the data type of the parameters in the collection
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param values
	 *            the values the names stand for, this array must be of the same
	 *            length as <code>keys</code>, a key at a certain position
	 *            belongs to the value at the same position in the array
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if the parameter is required,
	 *            <code>false</code> otherwise
	 * 
	 * @throws InconsistentCollectionException
	 *             if the lengths of <code>keys</code> and <code>values</code>
	 *             are different or the collection is inconsistent for some
	 *             other reason
	 * @throws IllegalValueException
	 *             if one of the values in <code>values</code> is not of type
	 *             <code>datatype</code>
	 * @throws DatatypeNotValidException
	 *             if the <code>datatype</code> is not one of the allowed values
	 */
	public MultiSelectionParameter(DataType datatype, String[] keys, Object[] values, String name, String comment, boolean required)
			throws InconsistentCollectionException, IllegalValueException, DatatypeNotValidException {
		super(datatype, keys, values, name, comment, required);
	}

	/**
	 * Constructor for a {@link MultiSelectionParameter}. The first
	 * option in the selection is selected by default.
	 * 
	 * @param datatype
	 *            the data type of the parameters in the collection
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param values
	 *            the values the names stand for, this array must be of the same
	 *            length as <code>keys</code>, a key at a certain position
	 *            belongs to the value at the same position in the array
	 * @param comments
	 *            the comments on the values in the collection
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if the parameter is required,
	 *            <code>false</code> otherwise
	 * 
	 * @throws InconsistentCollectionException
	 *             if the lengths of <code>keys</code> and <code>values</code>
	 *             are different or the collection is inconsistent for some
	 *             other reason
	 * @throws IllegalValueException
	 *             if one of the values in <code>values</code> is not of type
	 *             <code>datatype</code>
	 * @throws DatatypeNotValidException
	 *             if the <code>datatype</code> is not one of the allowed values
	 */
	public MultiSelectionParameter(DataType datatype, String[] keys, Object[] values, String[] comments, String name, String comment,
			boolean required) throws InconsistentCollectionException, IllegalValueException, DatatypeNotValidException {
		super(datatype, keys, values, comments, name, comment, required);
	}

	/**
	 * Creates a new {@link MultiSelectionParameter} from an array of
	 * {@link ParameterSet}s. The first option in the selection is selected by
	 * default.
	 * 
	 * @param values
	 *            the options/values in the collection
	 * @param name
	 *            the name of this {@link MultiSelectionParameter}
	 * @param comment
	 *            the comment on this {@link MultiSelectionParameter}
	 * @param required
	 *            <code>true</code> if this
	 *            {@link MultiSelectionParameter} is required,
	 *            <code>false</code> otherwise
	 */
	public MultiSelectionParameter(String name, String comment, boolean required, ParameterSet... values) {
		super( name, comment, required, values );
	}
	
	/**
	 * Creates a new {@link MultiSelectionParameter} from an array of
	 * {@link Class}es of {@link ParameterSet}s. The first option in the selection is selected by
	 * default. It is the lazy evaluation variant of {@link #MultiSelectionParameter(String, String, boolean, ParameterSet...)},
	 * i.e., the {@link ParameterSet}s are only created if necessary.
	 * 
	 * @param values
	 *            the array of {@link Class}es of {@link ParameterSet}s
	 * @param name
	 *            the name of this {@link MultiSelectionParameter}
	 * @param comment
	 *            the comment on this {@link MultiSelectionParameter}
	 * @param required
	 *            <code>true</code> if this
	 *            {@link MultiSelectionParameter} is required,
	 *            <code>false</code> otherwise
	 */
	public MultiSelectionParameter(String name, String comment, boolean required, Class<? extends ParameterSet>... values) {
		super( name, comment, required, values );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link MultiSelectionParameter} from its XML
	 * representation
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public MultiSelectionParameter(StringBuffer representation) throws NonParsableException {
		super(representation);
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.AbstractSelectionParameter#clone()
	 */
	@Override
	public MultiSelectionParameter clone() throws CloneNotSupportedException {
		MultiSelectionParameter clone = (MultiSelectionParameter) super.clone();
		clone.defaultSelected = defaultSelected.clone();
		clone.selected = selected.clone();
		return clone;

	}

	/**
	 * Sets the selection of the option with key <code>key</code> to the value
	 * of <code>selected</code>.
	 * 
	 * @param key
	 *            the key of the option
	 * @param selected
	 *            the selection value
	 * 
	 * @return <code>true</code> if the key could be found and set,
	 *         <code>false</code> otherwise
	 */
	public boolean setSelected(String key, boolean selected) {
		for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
			if (parameters.getParameterAt(i).getName().equals(key)) {
				this.selected[i] = selected;
				return true;
			}
		}
		return false;
	}

	/**
	 * Returns the selection value of the option with key <code>key</code>.
	 * 
	 * @param key
	 *            the key of the option
	 * 
	 * @return the selection value or <code>false</code> if no such option
	 *         exists
	 */
	public boolean isSelected(String key) {
		for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
			if (parameters.getParameterAt(i).getName().equals(key)) {
				return this.selected[i];
			}
		}
		return false;
	}

	/**
	 * Returns the value for the option with key <code>key</code>.
	 * 
	 * @param key
	 *            the key of the option
	 * 
	 * @return the value or <code>null</code> if the corresponding option either
	 *         does not exist or is not selected
	 */
	public Object getValueFor(String key) {
		for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
			if (parameters.getParameterAt(i).getName().equals(key)) {
				if (this.selected[i]) {
					return parameters.getParameterAt(i).getValue();
				} else {
					return null;
				}
			}
		}
		return null;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.AbstractSelectionParameter#checkValue(java.lang.Object)
	 */
	@Override
	public boolean checkValue(Object value) {
		if( value.getClass().isArray() ) {
			Object[] o = (Object[]) value;
			int i = 0;
			while( i < o.length && super.checkValue( o[i] ) );
			return i < o.length;
		} else {
			return super.checkValue( value );
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#setValue(java.lang.Object)
	 */
	@Override
	public void setValue(Object value) throws IllegalValueException {
		if( !value.getClass().isArray() ) {
			value = new Object[] { value };
		}
		Object[] o = (Object[]) value;
		int[] idx = new int[o.length];
		int i = 0;
		while( i < o.length && (idx[i] = check( o[i] )) >= 0 ) { 
			i++;
		}
		if( i < o.length ) {
			throw new IllegalValueException( name,"Problem with the " + i +"-th value: " + errorMessage );
		} else {
			Arrays.fill( selected, false );
			for( i = 0; i < idx.length; i++ ) {
				selected[idx[i]] = true;

				if( o[i] instanceof ParameterSetContainer ) {
					parameters.getParameterAt( idx[i] ).setValue( ((ParameterSetContainer)o[i]).getValue() );
				} else if( o[i] instanceof ParameterSet ) {
					parameters.getParameterAt( idx[i] ).setValue( o[i] );
				}
			}
			userSelected = true;
		}
	}

	/**
	 * Returns the indexes of the selected options.
	 * @return the indexes of the selected options.
	 */
	public int[] getSelected() {
		IntList li = new IntList();
		for(int i=0;i<selected.length;i++){
			if(selected[i]){
				li.add( i );
			}
		}
		return li.toArray();
	}

	/**
	 * Sets the selection of option with no. <code>idx</code> to
	 * <code>selected</code>.
	 * 
	 * @param idx
	 *            the index of the option
	 * @param selected
	 *            the selection value
	 * 
	 * @return <code>true</code> if the option exists and could be set,
	 *         <code>false</code> otherwise
	 */
	public boolean setSelected(int idx, boolean selected) {
		if (idx < this.selected.length) {
			this.selected[idx] = selected;
			return true;
		} else {
			return false;
		}
	}

	
	@Override
	public boolean isSelected(int idx) {
		return idx < selected.length && selected[idx];
	}

	/**
	 * Returns the value of the option with no. <code>idx</code>.
	 * 
	 * @param idx
	 *            the index of the option
	 * 
	 * @return the value or <code>null</code> if the corresponding option either
	 *         does not exists or is not selected
	 */
	public Object getValueFor(int idx) {
		if (idx < selected.length && selected[idx]) {
			return parameters.getParameterAt(idx).getValue();
		} else {
			return null;
		}
	}

	/**
	 * Returns the values of all selected options as an array.
	 * 
	 * @return the values of all selected options
	 */
	public Object[] getValues() {
		int count = 0;
		for (int i = 0; i < selected.length; i++) {
			if (selected[i]) {
				count++;
			}
		}
		Object[] values = new Object[count];
		count = 0;
		for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
			if (selected[i]) {
				values[count++] = parameters.getParameterAt(i).getValue();
			}
		}
		return values;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#getValue()
	 */
	@Override
	public Object getValue() {
		if (current > -1) {
			return parameters.getParameterAt(current).getValue();
		} else {
			return null;
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#hasDefaultOrIsSet()
	 */
	@Override
	public boolean hasDefaultOrIsSet() {
		int numSelected = 0;
		for (int i = 0; i < selected.length; i++) {
			if (selected[i]) {
				numSelected++;
				if (!parameters.getParameterAt(i).hasDefaultOrIsSet()) {
					errorMessage = "Parameter "+ parameters.getParameterAt(i).getName() + " not set.";
					return false;
				}
			}
		}
		if (numSelected < 1 && isRequired()) {
			return false;
		} else {
			return true;
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#reset()
	 */
	@Override
	public void reset() {
		selected = defaultSelected.clone();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.AbstractSelectionParameter#setDefault(java.lang.Object)
	 */
	@Override
	public void setDefault(Object defaultValue) throws IllegalValueException {
		if( selected == null ) {
			selected = new boolean[parameters.getNumberOfParameters()];
			defaultSelected = selected.clone();
		}
		setValue(defaultValue);
		System.arraycopy( selected, 0, defaultSelected, 0, selected.length );
		userSelected = false;
	}
	
	/**
	 * Sets the default selection of this {@link MultiSelectionParameter} to
	 * <code>defaultSelection</code>. This method also sets the current
	 * selection of this {@link MultiSelectionParameter} to the default
	 * and resets it such that {@link AbstractSelectionParameter#isUserSelected()}
	 * returns <code>false</code>.
	 * 
	 * @param defaultSelection
	 *            the new default selection
	 */
	public void setDefaultSelected(int[] defaultSelection){
		Arrays.fill( defaultSelected, false );
		for(int i=0;i<defaultSelection.length;i++){
			defaultSelected[defaultSelection[i]] = true;
		}
		System.arraycopy( defaultSelected, 0, selected, 0, selected.length );
		userSelected = false;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "MultiSelectionParameter";
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.AbstractSelectionParameter#appendFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer sup ) {
		super.appendFurtherInfos( sup );
		
		XMLParser.appendObjectWithTags(sup, selected, "selected");
		XMLParser.appendObjectWithTags(sup, defaultSelected, "defaultSelected");
		XMLParser.appendObjectWithTags(sup, current, "current");
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.AbstractSelectionParameter#extractFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos(StringBuffer representation) throws NonParsableException {
		super.extractFurtherInfos( representation );
		
		selected = XMLParser.extractObjectForTags(representation, "selected", boolean[].class );
		defaultSelected = XMLParser.extractObjectForTags(representation, "defaultSelected", boolean[].class );
		current = XMLParser.extractObjectForTags(representation, "current", int.class );
	}

	/**
	 * Returns the number of calls of
	 * {@link MultiSelectionParameter#next()} that can be called
	 * before <code>false</code> is returned.
	 * 
	 * @param afterIdx
	 *            the index after which shall be counted
	 * 
	 * @return the number of calls of
	 *         {@link MultiSelectionParameter#next()}
	 */
	public int getNumberOfNexts(int afterIdx) {
		int count = 0;
		for (int i = afterIdx + 1; i < selected.length; i++) {
			if (selected[i]) {
				count++;
			}
		}
		return count;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.AbstractCollectionParameter#hasDefault()
	 */
	@Override
	public boolean hasDefault() {
		int i = 0;
		while( i < defaultSelected.length && !defaultSelected[i] ) {
			i++;
		}
		return i < defaultSelected.length;
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		
		String selected = XMLParser.extractForTag( command, namePrefix ).toString();
		this.setValue( selected );
		if(this.getValue() instanceof GalaxyConvertible){
			//TODO change as soon as Galaxy supports multiple selections in ifs
			((GalaxyConvertible)this.getValue()).fromGalaxy( namePrefix+"_opt"+getSelected()[0], command );
		}
	}

	@Override
	protected boolean isDefaultSelection(int i) {
		return defaultSelected[i];
	}
}