package projects.encodedream.tools;
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
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Date;
import java.util.LinkedList;
import java.util.zip.GZIPOutputStream;

import de.jstacs.DataType;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import projects.encodedream.DeriveLabels;

public class DeriveLabelTool implements JstacsTool {

	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<>();
		
		pars.add(new FileParameter("Conservative peaks", "NarrowPeak file containing the conservative peaks", "narrowPeak,bed", true));
		pars.add(new FileParameter("Relaxed peaks", "NarrowPeak file containing the relaxed peaks", "narrowPeak,bed", true));
		pars.add(new FileParameter("FAI of genome", "FastA index file of the genome", "fai", true));
		try {
			pars.add(new SimpleParameter(DataType.INT, "Bin width", "The width of the genomic bins considered", true,new NumberValidator<Integer>(1, 10000),50));
			pars.add(new SimpleParameter(DataType.INT, "Region width", "The width of the genomic regions considered for overlaps", true,new NumberValidator<Integer>(1, 10000),50));
		} catch (ParameterException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		boolean bySummit = false;
		String consFile = ((FileParameter)parameters.getParameterAt(0)).getFileContents().getFilename();
		String relFile = ((FileParameter)parameters.getParameterAt(1)).getFileContents().getFilename();
		String faiFile = ((FileParameter)parameters.getParameterAt(2)).getFileContents().getFilename();
		int bin = (int) parameters.getParameterAt(3).getValue();
		int reg = (int) parameters.getParameterAt(4).getValue();
		
		File outfile = File.createTempFile("labels", ".temp.gz", new File("."));
		outfile.deleteOnExit();
		
		PrintWriter wr = new PrintWriter(new GZIPOutputStream(new FileOutputStream(outfile)));
		
		DeriveLabels.run(consFile, relFile, faiFile, wr, bin, reg, bySummit);//TODO 200bp bins???
		
		wr.close();
		TextResult tr = new TextResult("Labels", "Labels derived from peak files", new FileParameter.FileRepresentation(outfile.getAbsolutePath()), "tsv.gz", getToolName(), null, true);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
	}

	@Override
	public String getToolName() {
		return "Derive labels";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "labels";
	}

	@Override
	public String getDescription() {
		return "derives labels (S: summit, B: bound, A: ambiguous, U: unbound) for genomic regions based on ChIP-seq peak files";
	}

	@Override
	public String getHelpText() {
		return "**Derive labels** computes labels for genomic regions based on ChIP-seq peak files. The input ChIP-seq peak files must be provided in narrowPeak format"
				+ " and may come in *conservative*, i.e., IDR-thresholded, and *relaxed* flavors. In case only a single peak file is available, both of the corresponding parameters may be set "
				+ "to this one peak file. The parameter for the bin width defines the resolution of genomic regions that is assigned a label, while the parameter for the region width defines "
				+ "the size of the regions considered. If, for instance, the bin width is set to 50 and the region width to 100, regions of 100 bp shifted by 50 bp along the genome are labeled. "
				+ "The labels assigned may be *S* (summit) is the current bin contains the annotated summit of a conservative peak, *B* (bound) if the current region overlaps a conservative peak by at least "
				+ "half the region width, *A* (ambiguous) if the current region overlaps a relaxed peak by at least 1 bp, or *U* (unbound) if it overlaps with none of the peaks. The output is provided as a gzipped file "
				+ "*Labels.tsv.gz* with columns chromosome, start position, and label. This output file together with a protocol of the tool run is saved to the specified output directory.";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

}
