package projects.gemoma.JunitTest;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.Date;

import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.ui.cli.CLI.SysProtocol;
import projects.gemoma.ExtractRNAseqEvidence;
import projects.gemoma.ExtractRNAseqEvidence.Stranded;
import projects.gemoma.Extractor;
import projects.gemoma.Tools.Ambiguity;

/**
 * 
 * @author Jens Keilwagen
 */
public class GeMoMaJUnitTest {

//for creating infrastructure, clean up, and nice things
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		System.out.println("start: " + new Date());
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
		System.out.println("end: " + new Date());
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
		System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	}

//own test implementation	
	static String in = "projects/gemoma/JunitTest/given/";

	void compare( String name, String givenName, ResultSet[] rs ) throws IOException {
		TextResult tr = (TextResult) (rs[0].getResultForName(name));
		FileRepresentation fr = tr.getValue();
		assertTrue( name + ": The files differ!", FileUtils.contentEquals( new File(in + givenName), new File(fr.getFilename()) ));
	}
	
	@Test
	public void testExtractor() throws Exception {
		Extractor e = new Extractor(-1);
		
		//parameters
		ParameterSet ps = e.getToolParameters();
		ps.getParameterForName("annotation").setValue( in + "Arabidopsis_lyrata.v.1.0.31.chr.gff3");
		ps.getParameterForName("genome").setValue( in + "Arabidopsis_lyrata.v.1.0.31.dna.genome.fa");
		ps.getParameterForName("repair").setValue( true );
		ps.getParameterForName("full-length").setValue( false );
		ps.getParameterForName("Ambiguity").setValue( Ambiguity.AMBIGUOUS );
		ps.getParameterForName("proteins").setValue( true );
			
		//get results
		ResultSet[] rs = e.run(ps, new SysProtocol(), new ProgressUpdater(), 1).getValue();
		
		//compare that everything is correct
		assertTrue( "Differ number of ResultSets", rs.length == 1 );
		assertTrue( "Differ number of Results", rs[0].getNumberOfResults() == 3 );
		compare("cds-parts", "cds-parts.fasta", rs );
		compare("assignment", "assignment.tabular", rs );
		compare("proteins", "proteins.fasta", rs );
	}
	
	@Test
	public void testERE() throws Exception {
		ExtractRNAseqEvidence e = new ExtractRNAseqEvidence();
		
		//parameters
		ParameterSet ps = e.getToolParameters();
		ps.getParameterForName("Stranded").setValue( Stranded.FR_UNSTRANDED );
		ParameterSetContainer psc = (ParameterSetContainer) ps.getParameterAt(1);
		ExpandableParameterSet eps = (ExpandableParameterSet) psc.getValue();
		SimpleParameterSet sps = (SimpleParameterSet) ((ParameterSetContainer) eps.getParameterAt(0)).getValue();
		sps.getParameterAt(0).setValue(in + "RNAseq.bam");
		ps.getParameterForName("coverage output").setValue( true );

		//get results
		ResultSet[] rs = e.run(ps, new SysProtocol(), new ProgressUpdater(), 1).getValue();
		
		//compare that everything is correct
		assertTrue( "Differ number of ResultSets", rs.length == 1 );
		assertTrue( "Differ number of Results", rs[0].getNumberOfResults() == 2 );
		compare("coverage", "coverage.bedgraph", rs );
		compare("introns", "at-introns.gff", rs );
	}
}
