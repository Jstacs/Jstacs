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

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * A class for a {@link ParameterSet} that can be expanded by additional
 * {@link Parameter}s at runtime. The {@link Parameter}s are all of type
 * {@link ParameterSetContainer}. The contents of these {@link Parameter}s are
 * provided by a {@link ParameterSet}-template, which is cloned each time
 * {@link #addParameterToSet()} is called and set as value of a
 * {@link ParameterSetContainer} that is added to the {@link Parameter}s of the
 * {@link ExpandableParameterSet}. Already added {@link Parameter}s can also be
 * removed from the set using {@link #removeParameterFromSet()}.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class ExpandableParameterSet extends ParameterSet {

	/**
	 * The template for each {@link ParameterSet}
	 */
	protected ParameterSet template;
	/**
	 * A template for the name of the enclosing {@link ParameterSetContainer}
	 */
	protected String nameTemplate;
	/**
	 * A template for the comment of the enclosing {@link ParameterSetContainer}
	 */
	protected String commentTemplate;

	private int count;

	private int initCount, minCount, maxCount;

	/**
	 * Creates a new {@link ExpandableParameterSet} from a {@link Class} that
	 * can be instantiated using this {@link ExpandableParameterSet} and
	 * templates for the {@link ParameterSet} in each element of the array, the
	 * name and the comment that are displayed for the
	 * {@link ParameterSetContainer}s enclosing the {@link ParameterSet}s.
	 * 
	 * @param template
	 *            the template of the {@link ParameterSet}
	 * @param nameTemplate
	 *            the name-template
	 * @param commentTemplate
	 *            the comment-template
	 * @throws CloneNotSupportedException if the template could not be cloned
	 */
	public ExpandableParameterSet(ParameterSet template, String nameTemplate,
			String commentTemplate) throws CloneNotSupportedException {
		this( template, nameTemplate, commentTemplate, 1 );
	}

	/**
	 * Creates a new {@link ExpandableParameterSet} from a {@link Class} that
	 * can be instantiated using this {@link ExpandableParameterSet} and
	 * templates for the {@link ParameterSet} in each element of the array, the
	 * name and the comment that are displayed for the
	 * {@link ParameterSetContainer}s enclosing the {@link ParameterSet}s.
	 * 
	 * @param template
	 *            the template of the {@link ParameterSet}
	 * @param nameTemplate
	 *            the name-template
	 * @param commentTemplate
	 *            the comment-template
	 * @param initCount
	 *            the number of initial copies of the template
	 * @throws CloneNotSupportedException if the template could not be cloned
	 */
	public ExpandableParameterSet(ParameterSet template, String nameTemplate,
			String commentTemplate, int initCount) throws CloneNotSupportedException {
		this(template, nameTemplate, commentTemplate,initCount,initCount,Integer.MAX_VALUE);
	}
	
	/**
	 * Creates a new {@link ExpandableParameterSet} from a {@link Class} that
	 * can be instantiated using this {@link ExpandableParameterSet} and
	 * templates for the {@link ParameterSet} in each element of the array, the
	 * name and the comment that are displayed for the
	 * {@link ParameterSetContainer}s enclosing the {@link ParameterSet}s.
	 * 
	 * @param template
	 *            the template of the {@link ParameterSet}
	 * @param nameTemplate
	 *            the name-template
	 * @param commentTemplate
	 *            the comment-template
	 * @param initCount
	 *            the number of initial copies of the template
	 * @param minCount
	 *            the minimal number of copies of the template
	 * @param maxCount
	 *            the maximal number of copies of the template
	 * @throws CloneNotSupportedException if the template could not be cloned
	 */
	public ExpandableParameterSet(ParameterSet template, String nameTemplate,
			String commentTemplate, int initCount, int minCount, int maxCount ) throws CloneNotSupportedException {
		super();
		this.template = template;
		this.nameTemplate = nameTemplate;
		this.commentTemplate = commentTemplate;
		this.count = 0;
		this.initCount = initCount;
		for (int i = 0; i < initCount; i++) {
			addParameterToSet();
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ExpandableParameterSet} from its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public ExpandableParameterSet(StringBuffer representation)
			throws NonParsableException {
		super(representation);
	}

	/**
	 * Creates a new {@link ExpandableParameterSet} from a {@link ParameterSet}
	 * -array. The first element of this array is taken as a template for
	 * further calls of {@link #addParameterToSet()}.
	 * 
	 * @param templateAndContent
	 *            the content (and template)
	 * @param nameTemplate
	 *            the name-template
	 * @param commentTemplate
	 *            the comment-template
	 */
	public ExpandableParameterSet(ParameterSet[] templateAndContent,
			String nameTemplate, String commentTemplate) {
		if (templateAndContent.length == 0) {
			throw new IllegalArgumentException(
					"You must provide at least one ParameterSet.");
		}
		this.template = templateAndContent[0];
		this.nameTemplate = nameTemplate;
		this.commentTemplate = commentTemplate;

		if (notAllGivenParameterSetsAreOfTemplateType(templateAndContent)) {
			throw new IllegalArgumentException(
					"At least one of the given ParameterSets is not of "
							+ "the specified template-Type");
		}
		for (int i = 0; i < templateAndContent.length; i++) {
			parameters.add(new ParameterSetContainer(nameTemplate + " no. "
					+ (i + 1), commentTemplate, templateAndContent[i]));
		}
		this.count = templateAndContent.length;
		this.initCount = this.count;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ParameterSet#clone()
	 */
	@Override
	public ExpandableParameterSet clone() throws CloneNotSupportedException {
		ExpandableParameterSet clone = (ExpandableParameterSet) super.clone();
		clone.template = template.clone();
		return clone;
	}


	/**
	 * Adds a new {@link ParameterSetContainer} containing a clone of the
	 * {@link ParameterSet}-template to the set of {@link Parameter}s.
	 * 
	 * @throws CloneNotSupportedException
	 *             if the template could not be cloned
	 */
	public boolean addParameterToSet() throws CloneNotSupportedException {
		if( count < maxCount ) {
			ParameterSetContainer simplePar = new ParameterSetContainer(
					nameTemplate + " no. " + (++count), commentTemplate, template
							.clone());
			parameters.add(simplePar);
			return true;
		}
		return false;
	}

	/**
	 * First removes all previous added {@link ParameterSetContainer}s and
	 * afterwards adds all given {@link ParameterSet}s (in the given order)
	 * enclosed in new {@link ParameterSetContainer}s.
	 * 
	 * @param paramSetArray
	 *            the {@link ParameterSet}s to be set
	 * 
	 * @return <code>true</code>, if all given {@link ParameterSet}s are of the
	 *         same class as the defined template-{@link ParameterSet} of this
	 *         {@link ExpandableParameterSet}, <code>false</code> otherwise (In
	 *         this case, neither the previous added
	 *         {@link ParameterSetContainer}s are removed nor is any given
	 *         {@link ParameterSet} added. Hence the
	 *         {@link ExpandableParameterSet} stays unchanged.)
	 */
	public boolean replaceContentWith(ParameterSet[] paramSetArray) {

		if (notAllGivenParameterSetsAreOfTemplateType(paramSetArray)) {
			return false;
		}

		initParameterList(paramSetArray.length);
		for (int i = 0; i < paramSetArray.length; i++) {
			parameters.add(new ParameterSetContainer(nameTemplate + " no. "
					+ (i + 1), commentTemplate, paramSetArray[i]));
		}
		this.count = paramSetArray.length;

		return true;
	}

	private boolean notAllGivenParameterSetsAreOfTemplateType(
			ParameterSet[] temp) {

		Class clazz = template.getClass();

		for (int i = 0; i < temp.length; i++) {
			if (temp[i].getClass() != clazz)
				return true;
		}

		return false;
	}

	/**
	 * Returns <code>true</code> if there is still a {@link Parameter} that can
	 * be removed from the set.
	 * 
	 * @return
	 *         <code>true</code> if a {@link Parameter} can be removed from the set
	 */
	public boolean parameterRemovable() {
		return count > minCount;
	}

	/**
	 * Removes the last {@link Parameter} from set.
	 */
	public boolean removeParameterFromSet() {

		if (count > minCount) {
			parameters.remove(parameters.size() - 1);
			count--;
			return true;
		}
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ParameterSet#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags(buf, "superParameters");
		XMLParser.appendObjectWithTags(buf, template, "template");
		XMLParser.appendObjectWithTags(buf, nameTemplate, "nameTemplate");
		XMLParser.appendObjectWithTags(buf, commentTemplate, "commentTemplate");
		XMLParser.appendObjectWithTags(buf, initCount, "initCount");
		XMLParser.appendObjectWithTags(buf, minCount, "minCount");
		XMLParser.appendObjectWithTags(buf, maxCount, "maxCount");
		XMLParser.addTags(buf, "expandableParameterSet");
		return buf;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ParameterSet#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag( representation, "expandableParameterSet" );
		super.fromXML( XMLParser.extractForTag( representation, "superParameters") );
		template = XMLParser.extractObjectForTags( representation, "template", ParameterSet.class );
		nameTemplate = XMLParser.extractObjectForTags( representation, "nameTemplate", String.class );
		commentTemplate = XMLParser.extractObjectForTags( representation, "commentTemplate", String.class );
		initCount = XMLParser.extractObjectForTags( representation, "initCount", int.class );
		try {
			minCount = XMLParser.extractObjectForTags( representation, "minCount", int.class );
			maxCount = XMLParser.extractObjectForTags( representation, "maxCount", int.class );
		} catch( NonParsableException e ) {
			minCount=0;
			maxCount=Integer.MAX_VALUE;
		}
	}

	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean addLine, int indentation ) throws Exception {
		namePrefix = namePrefix+"_"+nameTemplate.replaceAll( "\\s", "_" );
		int nextIndentation = XMLParser.nextIndentation(indentation);
		
		for(int i=0;i<this.getNumberOfParameters()-count;i++){
			((GalaxyConvertible)getParameterAt( i )).toGalaxy( namePrefix+"_ps", configPrefix, depth+1, descBuffer, configBuffer, false, nextIndentation );
		}
		
		StringBuffer buf = new StringBuffer();
		StringBuffer buf2 = new StringBuffer();
		buf2.append( "$len($"+configPrefix+namePrefix+")" );
		XMLParser.addTags( buf2, namePrefix+"_len" );
		buf2.append( "\n" );
		buf2.append( "#for $"+namePrefix+"_i, $"+namePrefix+"_run in enumerate($"+configPrefix+namePrefix+")\n" );
		
		template.toGalaxy( namePrefix, namePrefix+"_run.", depth+1, buf, buf2, false, nextIndentation );
		
		buf2.append( "#end for" );
		
		XMLParser.addTagsAndAttributes( buf, "repeat", "name=\""+namePrefix+"\" title=\""+nameTemplate+"\" default=\""+initCount+"\" min=\""+minCount+"\" max=\""+maxCount+"\"", indentation );
		buf.insert(0, "\n");
		buf.append("\n");
		descBuffer.append( buf );
		configBuffer.append( buf2 );
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		namePrefix = namePrefix+"_"+nameTemplate.replaceAll( "\\s", "_" );
		while(count > 0){
			removeParameterFromSet();
		}
		int num = XMLParser.extractObjectForTags(  command, namePrefix+"_len", int.class );
		while(num > 0){
			addParameterToSet();
			num--;
		}
		int i=0;
		for(;i<getNumberOfParameters()-count;i++){
			((GalaxyConvertible)getParameterAt( i )).fromGalaxy( namePrefix+"_ps", command );
		}
		for(;i<getNumberOfParameters();i++){
			((GalaxyConvertible)getParameterAt( i )).fromGalaxy( namePrefix, command );
		}
	}

	public void toGalaxyTest( String namePrefix, int depth, StringBuffer testBuffer, int indentation ) throws Exception {
		namePrefix = namePrefix+"_"+nameTemplate.replaceAll( "\\s", "_" );
		int nextIndentation = XMLParser.nextIndentation(indentation);
		
		for( int i = 0; i < getNumberOfParameters(); i++ ) {
			StringBuffer buf = new StringBuffer();
			((GalaxyConvertible)getParameterAt(i)).toGalaxyTest( namePrefix, depth+1, buf, nextIndentation );
			if( buf.length()>0 ) {
				XMLParser.addTagsAndAttributes( buf, "repeat", "name=\""+namePrefix+"\"", indentation );
				testBuffer.append(buf);
			}
		}
	}	
}
