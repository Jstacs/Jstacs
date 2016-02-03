package de.jstacs.tools;

import de.jstacs.parameters.ParameterSet;
import de.jstacs.results.Result;

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

	public static class ResultEntry{
		
		private Class<? extends Result> clazz;
		private String format;
		private String name;
		
		public ResultEntry(Class<? extends Result> clazz, String format, String name) {
			super();
			this.clazz = clazz;
			this.format = format;
			this.name = name;
		}

		public Class<? extends Result> getDeclaredClass() {
			return clazz;
		}

		public String getFormat() {
			return format;
		}

		public String getName() {
			return name;
		}
		
		
		
		
	}
	
	/**
	 * Returns the input parameters of this tool. 
	 * The parameters should be empty but may have default values, which are used in all interface variants (command line, Galaxy, JavaFX GUI).
	 * @return the input parameters
	 */
	public ParameterSet getToolParameters();
	
	/**
	 * Runs the tool using the provided (now filled) parameters, which are in structure identical to those returned by {@link #getToolParameters()}. These parameters should only be used for this run and should not affect subsequent runs of the same tool.
	 * {@link Protocol} and {@link ProgressUpdater} may be used for indicating the tool's progress. Depending on the implementation and interface variant, these may be rendered differently (or not at all). This method returns all results of this tool encapsulated in a {@link ToolResult}.
	 * @param parameters the input parameters
	 * @param protocol the protocol
	 * @param progress the progress updater
	 * @return the results of this tool
	 * @throws Exception
	 */
	public ToolResult run(ParameterSet parameters, Protocol protocol, ProgressUpdater progress) throws Exception;

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
	 * May use <href="http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html">reStructureText</a> markup.
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
	 * Results are matched to default results by their name ({@link Result#getName()}) and result class, so the name supplied to {@link ResultEntry#ResultEntry(Class, String, String)} 
	 * must be identical to that of the final {@link Result} and should be unique (otherwise only one appearance of this name will be considered a default result) within a result class.
	 * The set of default results must always be returned in the same order.
	 * @return the default results (or <code>null</code> for no default results).
	 */
	public ResultEntry[] getDefaultResultInfos();
	
	//TODO getStandardResultInfos() //cf. https://wiki.galaxyproject.org/Admin/Tools/ToolConfigSyntax#A.3Cdata.3E_tag_set	
}
