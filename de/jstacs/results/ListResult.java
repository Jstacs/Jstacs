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

package de.jstacs.results;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;

import de.jstacs.AnnotatedEntity;
import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.ComparableElement;

/**
 * Class for a {@link Result} that contains a list or a matrix, respectively, of
 * {@link ResultSet}s. This class provides a way to build a hierarchy of
 * {@link Result}s and {@link ResultSet}s, or to create multi-dimensional
 * {@link Result}.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class ListResult extends Result {
	/**
	 * The internal list of {@link ResultSet}s that are part of this
	 * {@link ListResult}
	 */
	protected ResultSet[] list;

	private ResultSet annotation;
	
	private boolean export;


	/**
	 * Creates a new {@link ListResult} from an array of {@link ResultSet}s and
	 * a {@link ResultSet} of annotations, which may provide additional
	 * information on this {@link ListResult}.
	 * 
	 * @param name
	 *            the name of the {@link ListResult}
	 * @param comment
	 *            the comment on the {@link ListResult}
	 * @param annotation
	 *            an annotation on this {@link ListResult}
	 * @param results
	 *            the array of {@link ResultSet}s
	 */
	public ListResult( String name, String comment, ResultSet annotation, ResultSet... results ) {
		super(name, comment, DataType.LIST);
		this.list = new ResultSet[results.length];
		System.arraycopy(results, 0, list, 0, results.length);
		this.annotation = annotation;
		this.export = false;
	}
	
	/**
	 * Creates a new {@link ListResult} from a {@link java.util.Collection} of {@link ResultSet}s and
	 * a {@link ResultSet} of annotations, which may provide additional
	 * information on this {@link ListResult}.
	 * 
	 * @param name
	 *            the name of the {@link ListResult}
	 * @param comment
	 *            the comment on the {@link ListResult}
	 * @param annotation
	 *            an annotation on this {@link ListResult}
	 * @param coll
	 *            the {@link Collection} of {@link ResultSet}s
	 */
	public ListResult( String name, String comment, ResultSet annotation, Collection<ResultSet> coll ) {
		super(name, comment, DataType.LIST);
		this.list = new ResultSet[coll.size()];
		Iterator<ResultSet> it = coll.iterator();
		int i=0;
		while(it.hasNext()){
			list[i] = it.next();
			i++;
		}
		this.annotation = annotation;
		this.export = false;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ListResult} from the corresponding XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer}<code>representation</code> could
	 *             not be parsed
	 */
	public ListResult( StringBuffer representation ) throws NonParsableException {
		super(representation);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.results.Result#getResult()
	 */
	@Override
	public ResultSet[] getValue() {
		int i = 0, k;
		ResultSet[] resultsToShow = new ResultSet[list.length];
		for (; i < list.length; i++) {
			if (list[i] instanceof MeanResultSet) {
				ResultSet infos = ((MeanResultSet) list[i]).getInfos();
				NumericalResultSet statistics = ((MeanResultSet) list[i])
						.getStatistics();
				Result[] all = new Result[infos.getNumberOfResults()
						+ statistics.getNumberOfResults()];
				for (k = 0; k < infos.getNumberOfResults(); k++) {
					all[k] = infos.getResultAt(k);
				}
				for (k = 0; k < statistics.getNumberOfResults(); k++) {
					all[k + infos.getNumberOfResults()] = statistics
							.getResultAt(k);
				}
				resultsToShow[i] = new ResultSet(all);
			} else {
				resultsToShow[i] = list[i];
			}
		}
		return resultsToShow;
	}

	/**
	 * Returns a copy of the internal list of {@link ResultSet}s. The references
	 * to the {@link ResultSet}s in the array are not cloned.
	 * 
	 * @return the internal list of {@link ResultSet}s
	 */
	public ResultSet[] getRawResult() {
		ResultSet[] res = new ResultSet[list.length];
		System.arraycopy(list, 0, res, 0, res.length);
		return res;
	}

	/**
	 * Returns the number of {@link ResultSet}s in this {@link ListResult}
	 * @return the number of {@link ResultSet}s
	 */
	public int getNumberOfResultSets(){
		return list.length;
	}
	
	/**
	 * Returns a reference to the annotation of this {@link ListResult}.
	 * 
	 * @return the annotation of this {@link ListResult}
	 */
	public ResultSet getAnnotation() {
		return annotation;
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "listResult";
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#appendFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		if (annotation != null) {
			XMLParser.appendObjectWithTags(buf, annotation, "annotation");
		}
		XMLParser.appendObjectWithTags(buf, export, "export");
		XMLParser.appendObjectWithTags(buf, list, "list");
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#extractFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos( StringBuffer representation ) throws NonParsableException {
		try {
			annotation = XMLParser.extractObjectForTags( representation, "annotation", ResultSet.class );
		} catch (NonParsableException e) {
			annotation = null;
		}
		try{
			export = XMLParser.extractObjectForTags( representation, "export", Boolean.class );
		} catch (NonParsableException e){
			export = false;
		}
		list = (ResultSet[])XMLParser.extractObjectForTags(representation, "list" );
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringWriter s = new StringWriter();
		PrintWriter writer = new PrintWriter(s);
		print(writer);

		return s.toString();
	}

	/**
	 * Prints the information of this {@link ListResult} to the provided
	 * {@link PrintWriter}.
	 * 
	 * @param writer
	 *            the {@link PrintWriter}
	 */
	public void print(PrintWriter writer) {
		Result r;
		if (annotation != null) {
			DataType d;
			for (int i = 0; i < annotation.getNumberOfResults(); i++) {
				r = annotation.getResultAt(i);
				d = r.getDatatype();
				if (d != DataType.PNG && d != DataType.HTML
						&& d != DataType.LIST && d != DataType.STORABLE) {
					if (r.getName().equals("kind of assessment")) {
						writer.println("#");
					}
					writer.print("# ");
					writer.print(r.getName());
					writer.print(": ");
					writer.println(r.getValue().toString());
				}
			}
		}
		if (list != null) {
			//writer.println("# ");//TODO
			ResultSet[] res = getValue();
			boolean newNames;
			int i = 0, j, k;
			for (; i < res.length; i++) {
				newNames = i == 0;
				k = res[i].getNumberOfResults() - 1;
				// check if there are new names
				if (!newNames) {
					if (k + 1 != res[i - 1].getNumberOfResults()) {
						newNames = true;
					} else {
						for (j = 0; j <= k; j++) {
							if (!res[i].getResultAt(j).getName().equals(
									res[i - 1].getResultAt(j).getName())) {
								newNames = true;
								break;
							}
						}
					}
				}
				// only if we have new names
				if (newNames) {
					writer.print("# ");
					for (j = 0; j <= k; j++) {
						writer.print(res[i].getResultAt(j).getName());
						if (j < k) {
							writer.print("\t");
						} else {
							writer.println();
						}
					}

				}
				// write results
				for (j = 0; j <= k; j++) {
					writer.print(res[i].getResultAt(j).getValue());
					if (j < k) {
						writer.print("\t");
					} else {
						writer.println();
					}
				}
			}
			writer.flush();
		}
	}

	/**
	 * This method enables you to sort the entries of this container by a
	 * specified column.
	 * 
	 * @param columnName
	 *            the name of the column to be sorted
	 * 
	 * @return a new {@link ListResult}, where the entries of the specified
	 *         column are sorted
	 * 
	 * @throws IllegalArgumentException
	 *             if not all entries have a column with this name
	 */
	public ListResult sort(String columnName) throws IllegalArgumentException {
		ComparableElement[] c = new ComparableElement[list.length];
		Comparable comp = null;
		ResultSet r;
		int i = 0, k;
		for (; i < list.length; i++) {
			k = list[i].findColumn(columnName);
			if (k < 0) {
				if (list[i] instanceof MeanResultSet) {
					r = ((MeanResultSet) list[i]).getInfos();
					k = r.findColumn(columnName);
					if (k >= 0) {
						comp = (Comparable) r.getResultAt(k).getValue();
					}
				}

				if (k < 0) {
					throw new IllegalArgumentException(
							"Could not find such a column.");
				}
			} else {
				comp = (Comparable) list[i].getResultAt(k).getValue();
			}

			c[i] = new ComparableElement<ResultSet, Comparable>(list[i], comp);
		}
		Arrays.sort(c);
		ResultSet[] results = new ResultSet[list.length];
		for (i = 0; i < list.length; i++) {
			results[i] = (ResultSet) c[i].getElement();
		}
		return new ListResult(name, comment, annotation, results);
	}
	
	/**
	 * Returns if this {@link ListResult} is exported in {@link de.jstacs.tools.ui.galaxy.Galaxy}.
	 * @return if exported
	 */
	public boolean getExport() {
		return export;
	}

	/**
	 * Sets if this {@link ListResult} will be exported in {@link de.jstacs.tools.ui.galaxy.Galaxy}.
	 * @param export if exported
	 */
	public void setExport(boolean export) {
		this.export = export;
	}
	
	public boolean equalValues( AnnotatedEntity a ) {
		if( a instanceof ListResult ) {
			ListResult lr = (ListResult) a;
			if( list.length != lr.list.length ) {
				return false;
			} else {
				int i = 0;
				while( i < list.length && list[i].equals(lr.list[i]) ) {
					i++;				
				}
				return i == list.length;
			}
		} else {
			return false;
		}
	}
}
