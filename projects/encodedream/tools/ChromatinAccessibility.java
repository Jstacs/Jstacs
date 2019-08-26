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
import java.io.IOException;
import java.io.PrintStream;
import java.util.Date;
import java.util.LinkedList;
import java.util.zip.GZIPOutputStream;

import de.jstacs.DataType;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import projects.encodedream.Coverage;
import projects.encodedream.ObjectStream;
import projects.encodedream.Pileup;
import projects.encodedream.Pileup.CovPile;

public class ChromatinAccessibility implements JstacsTool {

	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<>();
		
		try {
			SelectionParameter sp = new SelectionParameter(DataType.PARAMETERSET, new String[]{"BAM/SAM","Bigwig"}, new ParameterSet[]{
					new SimpleParameterSet(new FileParameter("Input SAM/BAM", "The input file containing the mapped DNase-seq/ATAC-seq reads", "bam,sam", true)),
					new SimpleParameterSet(new FileParameter("Input Bigwig", "The input file containing the mapped DNase-seq/ATAC-seq reads", "bw,bigwig", true),new FileParameter("FastA index", "The genome index", "fai", true))
					
			}, "Data source", "The format of the input file containing the coverage information", true);
			
			pars.add(sp);
		} catch (Exception e1) {
			e1.printStackTrace();
		}

		try {
			pars.add(new SimpleParameter(DataType.INT, "Bin width", "The width of the genomic bins considered", true));
		} catch (DatatypeNotValidException e) {
			e.printStackTrace();
		}
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		File out = File.createTempFile("accessibility", ".temp.gz", new File("."));
		out.deleteOnExit();

		GZIPOutputStream os = new GZIPOutputStream(new FileOutputStream(out));


		SelectionParameter sp = (SelectionParameter) parameters.getParameterAt(0);
		int bin = (int) parameters.getParameterAt(1).getValue();

		if(sp.getSelected() == 0){
			String inputBAM = ((FileParameter)((ParameterSet)sp.getValue()).getParameterAt(0)).getFileContents().getFilename();	
			
			ObjectStream<CovPile> ps2 = new ObjectStream<>(10000);

			new Thread( ()->{
				try {
					Pileup.pileup(inputBAM, ps2,false,true);
					ps2.close();
				} catch (IOException | ArrayIndexOutOfBoundsException e) {
					e.printStackTrace();
					System.exit(1);
				}
			}).start();

			double lambdaBG = Coverage.estimateLambdaBG(ps2,inputBAM);
			//System.out.println(lambdaBG);

			ObjectStream<CovPile> ps = new ObjectStream<CovPile>(10000);
			new Thread( ()->{
				try {
					Pileup.pileup(inputBAM, ps, false, true);
					ps.close();
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(1);
				}
			}).start();


			Coverage.coverage(ps, lambdaBG, inputBAM, new PrintStream(os), bin);


		}else{
			String inputBW = ((FileParameter)((ParameterSet)sp.getValue()).getParameterAt(0)).getFileContents().getFilename();
			String faiFile = ((FileParameter)((ParameterSet)sp.getValue()).getParameterAt(1)).getFileContents().getFilename();
			
			Coverage.coverage(faiFile, inputBW, new PrintStream(os), bin);
		}
		
		os.close();
		
		TextResult tr = new TextResult("Chromatin accessibility", "Features computed from chromatin accessibility", new FileParameter.FileRepresentation(out.getAbsolutePath()), "tsv.gz", getToolName(), null, true);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
	}

	@Override
	public String getToolName() {
		return "Chromatin accessibility";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "access";
	}

	@Override
	public String getDescription() {
		return "computes chromatin accessibility features";
	}

	@Override
	public String getHelpText() {
		return "**Chromatin accessibility** computes several chromatin accessibility features from DNase-seq or ATAC-seq data provided as fold-enrichment tracks or SAM/BAM files of mapped reads. "
				+ "Features a computed with a certain resolution defined by the bin width parameter. Setting this parameter to 50, for instance, features are computed for non-overlapping 50 bp bins "
				+ "along the genome. If input data are provided as SAM/BAM file, coverage information is extracted and normalized locally in a similar fashion as proposed for the MACS peak caller. "
				+ "Output is provided as a gzipped file *Chromatin_accessibility.tsv.gz* with columns chromosome, start position of the bin, minimum coverage and median coverage in the current bin, "
				+ "minimum coverage in 1000 bp regions before and after the current bin, maximum coverage in 1000 bp regions before and after the current bin, the number of steps in the coverage profile, "
				+ "and the number of monotonically increasing and decreasing steps in the coverage profile of the current bin. "
				+ "This output file together with a protocol of the tool run is saved to the specified output directory.";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

	@Override
	public ToolResult[] getTestCases(String path) {
		return null;
	}

	@Override
	public void clear() {		
	}

	@Override
	public String[] getReferences() {
		return null;
	}
	
	
}
