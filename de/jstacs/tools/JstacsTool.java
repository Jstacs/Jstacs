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

package de.jstacs.tools;

import de.jstacs.DataType;
import de.jstacs.parameters.AbstractSelectionParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.results.Result;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.tools.ui.galaxy.Galaxy;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor;

/**
 * Interface for a generic Jstacs tool.
 * 
 * Implementations of this interface may be used in the {@link CLI}, {@link Galaxy}, and JavaFX-based GUI classes
 * for creating generic interfaces to tools created using Jstacs
 * 
 * @author Jan Grau
 *
 */
public interface JstacsTool {

	/**
	 * Representation of default results.
	 * 
	 * @see JstacsTool#getDefaultResultInfos()
	 * 
	 * @author Jan Grau
	 *
	 */
	public static class ResultEntry{
		
		private Class<? extends Result> clazz;
		private String format;
		private String name;
		
		/**
		 * Creates a new default result entry.
		 * @param clazz the class of the default result
		 * @param format the format of the default result; should conform to Galaxy formats; if <code>null</code>, format is inferred from <code>clazz</code>
		 * @param name the name of the result, used as identifier
		 * @see GalaxyAdaptor#getDefaultExtension(Class)
		 */
		public ResultEntry(Class<? extends Result> clazz, String format, String name) {
			super();
			this.clazz = clazz;
			this.format = format;
			this.name = name;
		}

		/**
		 * Returns the class declared for the default result.
		 * @return the class
		 */
		public Class<? extends Result> getDeclaredClass() {
			return clazz;
		}

		/**
		 * Returns the format of the result.
		 * @return the format
		 */
		public String getFormat() {
			return format;
		}

		/**
		 * Returns the name of the result.
		 * @return the name
		 */
		public String getName() {
			return name;
		}	
	}
	
	/**
	 * This method returns a short {@link String} representation of simple parameters.
	 * 
	 * @return a short {@link String} representation of simple parameters
	 */
	public static String getSimpleParameterInfo( ParameterSet parameters ) {
		String res = null;
		ParameterSet inner;
		for( int i = 0; i < parameters.getNumberOfParameters(); i++ ) {
			Parameter p = parameters.getParameterAt(i);
			inner = null;
			if( (p instanceof SimpleParameter || p instanceof AbstractSelectionParameter) ){ 
				if( res == null ) {
					res = "";
				} else {
					res +="; ";
				}
				Object o;
				if ( p instanceof SelectionParameter ) {
					SelectionParameter a = (SelectionParameter) p;
					o = a.getParametersInCollection().getParameterAt(a.getSelected()).getName();
				} else {
					o = p.getValue();
				}
				res += p.getName() + ": " + o;
				if ( p instanceof SelectionParameter && p.getDatatype() == DataType.PARAMETERSET ) {
					inner = (ParameterSet) p.getValue();
				}
			}
			if( p instanceof ParameterSetContainer ) {
				inner = ((ParameterSetContainer)p).getValue();
			}
			if( inner != null ) {
				String subRes = getSimpleParameterInfo( inner );
				if( subRes != null ) {
					if( res == null ) {
						res = "";
					} else {
						res +="; ";
					}
					res += subRes;
				}				
			} 
		}
		return res;
	}
	
	/**
	 * Returns the input parameters of this tool. 
	 * The parameters should be empty but may have default values, which are used in all interface variants (command line, Galaxy, JavaFX GUI).
	 * @return the input parameters
	 */
	public ToolParameterSet getToolParameters();
	
	/**
	 * Runs the tool using the provided (now filled) parameters, which are in structure identical to those returned by {@link #getToolParameters()}. These parameters should only be used for this run and should not affect subsequent runs of the same tool.
	 * {@link Protocol} and {@link ProgressUpdater} may be used for indicating the tool's progress. Depending on the implementation and interface variant, these may be rendered differently (or not at all). This method returns all results of this tool encapsulated in a {@link ToolResult}.
	 * @param parameters the input parameters
	 * @param protocol the protocol
	 * @param progress the progress updater
	 * @param threads the maximum number of threads that may be used for this run of the tool
	 * @return the results of this tool
	 * @throws Exception if the tool can not be run properly
	 */
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception;

	/**
	 * Returns a descriptive, human readable name for this tool.
	 * @return the name
	 */
	public String getToolName();
	
	/**
	 * Returns a descriptive, human readable version for this tool.
	 * @return the version
	 */
	public String getToolVersion();
	
	/**
	 * Returns a name (preferably short and without spaces) for referring to this tool on the command line.
	 * @return the short name
	 */
	public String getShortName();
	
	/**
	 * Returns a short description (half a sentence) on what this tool does.
	 * @return the description
	 */
	public String getDescription();
	
	/**
	 * Returns a detailed help text for this tool, describing the purpose of the tool, all parameters and results.
	 * May use <a href="http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html">reStructuredText</a> markup.
	 * @return the help text
	 */
	public String getHelpText();
	
	/**
	 * Returns {@link ResultEntry}s for the default results of this {@link JstacsTool}. These results must be a subset of
	 * all results that are returned by successful runs of this tool regardless of the input parameters.
	 * If no results should be defined as defaults, this method may return <code>null</code>.
	 * 
	 * Currently, this information is only used in the {@link de.jstacs.tools.ui.galaxy.Galaxy} environment, but may be used
	 * in the JavaFX GUI in the future for building workflows as well.
	 * 
	 * For each entry either the format or the result class should be supplied. If a format is specified this overrides the default formats
	 * for this {@link Result} type. Default formats for a given {@link Result} class may be obtained from 
	 * {@link de.jstacs.tools.ui.galaxy.GalaxyAdaptor#getDefaultExtension(Class)}. If this method returns <code>null</code>, the output format
	 * may be guessed by the environment (e.g., Galaxy).
	 * Results are matched to default results by their name ({@link Result#getName()}) and result class, so the name supplied to {@link JstacsTool.ResultEntry#JstacsTool.ResultEntry(Class, String, String)} 
	 * must be identical to that of the final {@link Result} and should be unique (otherwise only one appearance of this name will be considered a default result) within a result class.
	 * The set of default results must always be returned in the same order.
	 * @return the default results (or <code>null</code> for no default results).
	 */
	public ResultEntry[] getDefaultResultInfos();
	
}
