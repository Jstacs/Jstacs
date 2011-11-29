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

package de.jstacs.data.sequences.annotation;


/**
 * This {@link SequenceAnnotationParser} returns always <code>null</code> as {@link SequenceAnnotation}.
 * That is it add no {@link SequenceAnnotation} to each {@link de.jstacs.data.Sequence}.
 * 
 * <br><br>
 * 
 * The class should be used by accessing {@link NullSequenceAnnotationParser#DEFAULT_INSTANCE}.
 * 
 * @author Jens Keilwagen
 */
public final class NullSequenceAnnotationParser implements SequenceAnnotationParser {

	/**
	 * The only instance of this class which is publicly available. 
	 */
	public static final NullSequenceAnnotationParser DEFAULT_INSTANCE = new NullSequenceAnnotationParser();
	
	private NullSequenceAnnotationParser(){};
	
	public void addToAnnotation( String unparsed ) {
	}

	public void clearAnnotation() {
	}

	public SequenceAnnotation[] getCurrentAnnotation() {
		return null;
	}

	public String parseAnnotationToComment( char commentChar, SequenceAnnotation... annotations ) {
		return "" + commentChar;
	}
}
