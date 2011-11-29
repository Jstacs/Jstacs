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


/**
 * This is a simple class that extracts {@link String}s. It enables the user to
 * create a {@link de.jstacs.data.Sample} from a bunch of {@link String}s. It
 * might be useful if one likes to select or create the {@link String}s by an
 * own procedure. The extracted Strings which should be parsed to
 * {@link de.jstacs.data.Sequence} will not be annotated, since
 * {@link SimpleStringExtractor#getCurrentSequenceAnnotations()} always returns
 * <code>null</code>. 
 * 
 * @author Jens Keilwagen
 */
public class SimpleStringExtractor extends AbstractStringExtractor {

	private String[] strings;

	private int idx;

	/**
	 * This constructor packs the {@link String}s in an instance of
	 * {@link SimpleStringExtractor}.
	 * 
	 * @param strings
	 *            the {@link String}s that will be packed
	 * 
	 * @see AbstractStringExtractor#AbstractStringExtractor(char)
	 */
	public SimpleStringExtractor( String... strings ) {
		super( FASTA );
		this.strings = strings.clone();
		idx = 0;
	}

	/* (non-Javadoc)
	 * @see java.util.Enumeration#hasMoreElements()
	 */
	public boolean hasMoreElements() {
		return idx < strings.length;
	}

	/* (non-Javadoc)
	 * @see java.util.Enumeration#nextElement()
	 */
	public String nextElement() {
		return strings[idx++];
	}
}
