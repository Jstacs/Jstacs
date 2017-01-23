package projects.gemoma.JunitTest;

import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Date;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.ui.cli.CLI.SysProtocol;
import projects.gemoma.ExtractRNAseqEvidence;
import projects.gemoma.ExtractRNAseqEvidence.Stranded;
import projects.gemoma.Extractor;
import projects.gemoma.GeMoMa;
import projects.gemoma.GeMoMaAnnotationFilter;
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

	void assertFile( String name, String givenName, ResultSet rs ) throws AssertionError, IOException {
		long l = -100;
		String l1, l2;
		File f1, f2;
		BufferedReader r1, r2;
		l1=l2=null;
		f1=f2=null;
		r1=r2=null;
		try { 
			TextResult tr = (TextResult) (rs.getResultForName(name));
			FileRepresentation fr = tr.getValue();

			//File help = File.createTempFile(name.replaceAll(" ", "_"), ".tmp", new File(in));
			//FileManager.writeFile(help, fr.getContent() );
			f1 = new File(in + givenName);
			f2 = new File(fr.getFilename()); //help
			
			if( f1.exists() != f2.exists() ) {
				throw new AssertionError(name + ": The files differ.\nFile 1: "  + f1 );
			} else {
				r1= new BufferedReader( new FileReader(f1) );
				r2= new BufferedReader( new FileReader(f2) );
				l=0;
				while( true ) {
					l1=r1.readLine();
					l2=r2.readLine();
					
					if( l1==null || l2 == null || !l1.equals(l2) ) {
						break;
					}
					l++;
				}
				r1.close();
				r2.close();
				if( l1 == null && l2 == null ) {
					//okay
				} else {
					throw new RuntimeException();
				}
			}
			//help.deleteOnExit();
		} catch( Exception e ) {
			if( l < 0 ) {
				throw new AssertionError(name + ": The files differ."
						+ "\nFile 1 ("  + f1.getAbsolutePath() + "): " + (f1.exists() ? "exists" : "does not exist")
						+ "\nFile 1 ("  + f2.getAbsolutePath() + "): " + (f2.exists() ? "exists" : "does not exist"),
						e );
			} else {
				throw new AssertionError(name + ": The files differ in line "+ l+"."
						+ "\nFile 1("+f1.getAbsolutePath()+"): " + (l1==null?"[EOF]":l1) 
						+ "\nFile 2("+f2.getAbsolutePath()+"): " + (l2==null?"[EOF]":l2) );
			}
		} finally {
			if( r1 != null ) {
				r1.close();
			}
			if( r2 != null ) {
				r2.close();
			}
		}
		
		//assertTrue( name + ": The files differ!", FileUtils.contentEquals( f1, f2 ));
	}
	
	void assertStructure( ResultSet[] rs, int num, int... len ) {
		assertTrue( "Different number of ResultSets", rs.length == num );
		for( int i = 0; i < len.length; i++ ) {
			assertTrue( "Different number of Results in ResultSet " + i, rs[i].getNumberOfResults() == len[i] );
		}
	}
	
	void set( ParameterSet ps, int idx, String value ) throws IllegalValueException {
		ParameterSetContainer psc = (ParameterSetContainer) ps.getParameterAt(idx);
		ExpandableParameterSet eps = (ExpandableParameterSet) psc.getValue();
		SimpleParameterSet sps = (SimpleParameterSet) ((ParameterSetContainer) eps.getParameterAt(0)).getValue();
		sps.getParameterAt(0).setValue(in + value);
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
		assertStructure( rs, 1, 3 );
		assertFile("cds-parts", "cds-parts.fasta", rs[0] );
		assertFile("assignment", "assignment.tabular", rs[0] );
		assertFile("proteins", "proteins.fasta", rs[0] );
	}
	
	@Test
	public void testERE() throws Exception {
		ExtractRNAseqEvidence e = new ExtractRNAseqEvidence();
		
		//parameters
		ParameterSet ps = e.getToolParameters();
		ps.getParameterForName("Stranded").setValue( Stranded.FR_UNSTRANDED );
		set( ps, 1, "RNAseq.bam");
		ps.getParameterForName("coverage output").setValue( true );

		//get results
		ResultSet[] rs = e.run(ps, new SysProtocol(), new ProgressUpdater(), 1).getValue();
		
		//compare that everything is correct
		assertStructure( rs, 1, 2 );
		assertFile("coverage", "coverage.bedgraph", rs[0] );
		assertFile("introns", "at-introns.gff", rs[0] );
	}
	
	@Test
	public void testGeMoMa() throws Exception {
		GeMoMa e = new GeMoMa(-1, 3600, 60*60*24*7);
		
		//parameters
		ParameterSet ps = e.getToolParameters();
		ps.getParameterForName("tblastn results").setValue( in + "tblastn.tabular");
		ps.getParameterForName("contig threshold").setValue( 0.4 );
		ps.getParameterForName("predictions").setValue(10);
		ps.getParameterForName("target genome").setValue( in + "TAIR10_chr_all.fas");
		ps.getParameterForName("assignment").setValue( in + "assignment.tabular");
		ps.getParameterForName("cds parts").setValue( in + "cds-parts.fasta");
		ps.getParameterForName("query proteins").setValue( in + "proteins.fasta");
		ps.getParameterForName("introns").setValue( in + "at-introns.gff");
		SelectionParameter sp = (SelectionParameter) ps.getParameterForName("coverage");
		sp.setValue("UNSTRANDED");
		((SimpleParameterSet) sp.getValue()).getParameterAt(0).setValue( in + "coverage.bedgraph");
		
		
		//get results
		ResultSet[] rs = e.run(ps, new SysProtocol(), new ProgressUpdater(), 1).getValue();
		
		//compare that everything is correct
		assertStructure( rs, 1, 2 );
		assertFile("predicted annotation", "prediction.gff", rs[0] );
		assertFile("predicted protein", "predicted-protein.fasta", rs[0] );
	}
	
	@Test
	public void testGAF() throws Exception {
		GeMoMaAnnotationFilter e = new GeMoMaAnnotationFilter();
		
		//parameters
		ParameterSet ps = e.getToolParameters();
		set( ps, 6, "prediction.gff");
		
		//get results
		ResultSet[] rs = e.run(ps, new SysProtocol(), new ProgressUpdater(), 1).getValue();
		
		//compare that everything is correct
		assertStructure( rs, 1, 1 );
		assertFile("filtered predictions", "filtered_predictions.gff", rs[0] );
	}
}
