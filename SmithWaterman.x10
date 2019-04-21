package smithwatermanalgo;
import x10.lang.String;
import x10.io.Console;
import x10.io.File;
import x10.util.StringBuilder;
import x10.io.FileReader;
import x10.lang.Exception;
import x10.util.HashMap;
import x10.array.Array_2;


public class SmityWaterman {
    private val seq1:String;
    private val seq2:String;
    private val length1:Int;
    private val length2:Int;
    private var score:Array_2[Double];
    private val MATCH_SCORE:Double;
    private val GAPOPEN_SCORE:Double;
    private val GAPEXT_SCORE:Double;
    private val DR_LEFT:Int = 1;// 0001
    private val DR_UP:Int = 2;// 0010
    private val DR_DIAG:Int = 4;// 0100
    private val DR_ZERO:Int = 8;// 1000
    private val BLOSUM62:Array_2[Long];
    private val MAP:HashMap[String, Int];
    private var sizesOfVerticalGaps:Array[Short];
    private var sizesOfHorizontalGaps:Array[Short];
    private var prevCells:Array_2[Double];
    private val blosumReader: BlosumReader;
    

    public def this(fasta1:Rail[String], fasta2:Rail[String]) {
        this.seq1 = FastaReader.readFastaFile(fasta1);
        this.seq2 = FastaReader.readFastaFile(fasta2);
        this.blosumReader = new BlosumReader('./BLOSUM62');

    }


}

class class FastaReader {
    public static def readFastaFile(fastaFileName:String):String {
        val fastaFile:File = new File(fastaFileName);
        val fastaReader:FileReader = fastaFile.openRead();
        val header:String = fastaReader.readLine();
        val builder:StringBuilder = new StringBuilder();
        var line:String = null;

        while (true) {
            try {
                line = fastaReader.readLine().trim();
            } catch (e:EOFException) {
                break;
            }
            builder.add(line);
        }
        return builder.toString();
    } 
}

 class BlosumReader {
    private var BLOSUM62: Array_2[Int];
    private val SeqToNum: HashMap[String, Int];
    private val NumToSeq: HashMap[Int, String];
    private val NUMOFSEQ: Int = 23n;

    public def this(blosumFileName: String) {
        this.SeqToNum = this.initSeqToNum();
        this.NumToSeq = this.initNumToSeq();
        this.BLOSUM62 = new Array_2[Int](NUMOFSEQ + 1n, NUMOFSEQ + 1n);
        this.readBlosumFile(blosumFileName);
    }

    private def readBlosumFile(blosumFileName: String) {
        val fastaFile: File = new File(blosumFileName);
        val fastaReader: FileReader = fastaFile.openRead();
        val header: String = fastaReader.readLine().trim();
        var line: String = null;
        var chars: Rail[String] = new Rail[String](NUMOFSEQ + 1n);
        for (i in 0n .. (NUMOFSEQ)) {
            line = fastaReader.readLine().trim();
            chars = line.split(" ");
            for (j in 0n .. NUMOFSEQ) {
                this.BLOSUM62(i, j) = Int.parseInt(chars(j + 1n));
            }
        }
    }

    private def initSeqToNum(): HashMap[String, Int] {
        val map:HashMap[String, Int] = new HashMap[String, Int]();
        map.put("A", 0n);
        map.put("R", 1n);
        map.put("N", 2n);
        map.put("D", 3n);
        map.put("C", 4n);
        map.put("Q", 5n);
        map.put("E", 6n);
        map.put("G", 7n);
        map.put("H", 8n);
        map.put("I", 9n);
        map.put("L", 10n);
        map.put("K", 11n);
        map.put("M", 12n);
        map.put("F", 13n);
        map.put("P", 14n);
        map.put("S", 15n);
        map.put("T", 16n);
        map.put("W", 17n);
        map.put("Y", 18n);
        map.put("V", 19n);
        map.put("B", 20n);
        map.put("Z", 21n);
        map.put("X", 22n);
        return map;
    }

    private def initNumToSeq(): HashMap[Int, String] {
        val map:HashMap[Int, String] = new HashMap[Int, String]();
        map.put(0n, "A");
        map.put(1n, "R");
        map.put(2n, "N");
        map.put(3n, "D");
        map.put(4n, "C");
        map.put(5n, "Q");
        map.put(6n, "E");
        map.put(7n, "G");
        map.put(8n, "H");
        map.put(9n, "I");
        map.put(10n, "L");
        map.put(11n, "K");
        map.put(12n, "M");
        map.put(13n, "F");
        map.put(14n, "P");
        map.put(15n, "S");
        map.put(16n, "T");
        map.put(17n, "W");
        map.put(18n, "Y ");
        map.put(19n, "V");
        map.put(20n, "B");
        map.put(21n, "Z");
        map.put(22n, "X");
        return map;
    }

    public def getSeqToNum(): HashMap[String, Int] {
        return this.SeqToNum;
    }

    public def getNumToSeq(): HashMap[Int, String] {
        return this.NumToSeq;
    }

    public def getBlosum62(): Array_2[Int] {
        return this.BLOSUM62;
    }
}






