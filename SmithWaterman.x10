package smithwatermanalgo;
import x10.lang.String;
import x10.io.Console;
import x10.io.File;
import x10.util.StringBuilder;
import x10.io.FileReader;
import x10.lang.Exception;
import x10.util.HashMap;
import x10.array.Array_2;
import x10.lang.Math;
import x10.io.Console;


public class SmithWaterman {

    //two sequences of AA
    private val seq1:String;
    private val seq2:String;
    private val length1:Int;
    private val length2:Int;

    //the score matrix
    private var score:Array_2[Int];

    private var scoreLeft:Array_2[Int];
    private var scoreUp:Array_2[Int];

    //similarity fuction constant
    private val SCORE_THRESHOLD:Int = 19.9;
    private val MATCH_SCORE:Int;
    private val GAP_OPENING_PANALTY:Int;
    private val GAP_EXTENSION_PANALTY:Int;
    private val INDEL_SCORE:Int = -9;

    //constant in direction matrix
    private val DR_LEFT:Int = 1n;// 0001
    private val DR_UP:Int = 2n;// 0010
    private val DR_DIAG:Int = 4n;// 0100
    private val DR_ZERO:Int = 8n;// 1000



    //direction matrix
    private var prevCells:Array_2[Int];

    private val blosum62:Array_2[Int];

    private val seqToNum:HashMap[String, Int];
    private val blosumReader:BlosumReader;
    private val fastaReader:FastaReader;
    

    public def this(fastaName1:String, fastaName2:String, blosumFileName:String, openPanalty:Int, extensionPanalty:Int) {
        
        fastaReader = new FastaReader();
        this.blosumReader = new BlosumReader(blosumFileName);
        this.blosum62 = blosumReader.getBlosum62();
        this.seqToNum = blosumReader.getSeqToNum();

        this.seq1 = fastaReader.readFastaFile(fastaName1);
        this.seq2 = fastaReader.readFastaFile(fastaName2);
        this.length1 = seq1.length();
        this.length2 = seq2.length();

        this.GAP_OPENING_PANALTY = openPanalty;
        this.GAP_EXTENSION_PANALTY = extensionPanalty;

        score = new Array_2[Int](length1 + 1, length2 + 1);
        scoreLeft = new Array_2[Int](length1 + 1, length2 + 1);
        scoreUp = new Array_2[Int](length1 + 1, length2 + 1);
        prevCells = new Array_2[Int](length1 + 1, length2 + 1);

        buildMatrix();
    }

    private def similarity(i:Int, j:Int):Long {
        return blosum62(seqToNum.get(seq1.charAt(i-1n)), seqToNum.get(seq2.charAt(j-1n)));
    }

    public def buildMatrix() {


        //base case
        score(0, 0) = 0n;
        scoreLeft(0, 0) = 0n;
        scoreUp(0, 0) = 0n;
        prevCells(0,0) = DR_ZERO;

        //the first row
        for(i in 1..length1) {
            score(i, 0) = 0n;
            scoreLeft(i, 0) = 0n;
            scoreUp(i, 0) = 0n;
            prevCells(i, 0) = DR_ZERO;
        }

        //the first column
        for(j in 1..length2) {
            score(0, j) = 0n;
            scoreLeft(0, j) = 0n;
            scoreUp(0, j) = 0n;
            prevCells(0, j) = DR_ZERO;
        }

        diagnalCover(1n, 1n, length1, length2);

    }

    public def diagnalCover(var a1:Int, var b1:Int, var a2:Int, var b2:Int) {

        //start from row 0
        for(j in b1 .. b2) {
            var i:Int = a1;
            while(i <= a2 && j >= b1) {
                calculateScore(i, j);
                i++;
                j--;
            }
        }

        //continue from final col
        for(i in a1+1 .. a2) {
            var j:Int = b2;
            while(i <=a2 && j >= b1) {
                calculateScore(i, j);
                i++;
                j--;
            }
        }
    }


    public def calculateScore(var i:Int, var j:Int) {
        var diagScore:Int = score(i-1, j-1) + similarity(i, j);

        var newOpenGapLeftScore:Int = score(i, j-1) - GAP_OPENING_PANALTY;
        var newExtentionGapLeftScore:Int = scoreLeft(i, j-1) - GAP_EXTENSION_PANALTY;
        scoreLeft(i, j) = Math.max(newOpenGapLeftScore, newExtentionGapLeftScore);

        var newOpenGapUpScore:Int = score(i-1, j) - GAP_OPENING_PANALTY;
        var newExtentionGapUpScore:Int = scoreUp(i-1, j) - GAP_EXTENSION_PANALTY;
        scoreUp(i, j) = Math.max(newOpenGapUpScore,newExtentionGapUpScore);


        var upScore:Int = scoreUp(i, j-1);
        var leftScore:Int = scoreLeft(i-1, j);

        score(i, j) = Math.max(diagScore, Math.max(upScore, Math.max(leftScore, 0)));
        prevCells(i, j) = 0;

        if (diagScore == score(i, j)) {
            prevCells(i, j) |= DR_DIAG;
        }
        if (leftScore == score(i, j)) {
            prevCells(i, j) |= DR_LEFT;
        }
        if (upScore == score(i, j)) {
            prevCells(i, j) |= DR_UP;
        }
        if (0 == score(i, j)) {
            prevCells(i, j) |= DR_ZERO;
        }
    }

    public def getMaxScore():Int {
        var maxScore:Int = 0;

        for(i in 1 .. length1) {
            for(j in 1 .. length2) {
                if(score(i, j) > maxScore) {
                    maxScore = score(i, j);
                }
            }
        }
        return maxScore;
    }


    // TODO: printAlignments()

    //returns the end point of tracing back (the top left cell) and the number of matches
    public def traceback(var i:Int, var j:Int):Rail[Int] {
        var num:Int = 0;

        //find the direction to traceback
        while (true)
        {
            if ((prevCells(i, j) & DR_LEFT) > 0) {
                num ++;
                if (score(i-1, j)>0) i--;
                else    break;              
            }
            if ((prevCells(i, j) & DR_UP) > 0) {
                num ++;
//          return traceback(i, j-1);
                if (score(i, j-1)>0) j--;
                else    break;              
            }
            if ((prevCells(i, j) & DR_DIAG) > 0) {
                num ++;
//          return traceback(i-1, j-1);
                if (score(i-1, j-1)>0) {i--;j--;}
                else     break;             
            }
        }
        var point:Rail[Int] = [i, j, num];
        return point;
    }

    public def getMatchNumber():Int {

        var matches:Int = 0;

        for(i in 1..length1) {
            for(j in 1..length2) {
                if(score(i, j) > SCORE_THRESHOLD && score(i, j) > score(i-1, j)
                    && score(i, j) > score(i ,j-1) && score(i, j) > score(i-1, j-1))
                {
                    if(i == length1 || j == length2 || score(i, j) > score(i+1, j+1))
                    {
                        var endPoint:Rail[Int] = traceback(i, j);

                        matches += endPoint(2);
                    }
                }
            }
        }
    }

    public def getTototalNumber():Int {

    }



    def main(argv: Array[String](1)) {
        x10.Console.OUT.println("Input the FASTA_FILE_1 FASTA_FILE_2 MATCH_FILE GAP_OPENING_PANALTY GAP_EXTENSION_PANALTY");
        val s = x10.io.Console.IN.readLine();

        val param = s.split(" ");

        val fasta1:String = param(0);
        val fasta2:String = param(1);
        val match:String = param(2);
        val openPanalty:Int = param(3) as Int;
        val extPanalty:Int = param(4) as Int;

        val sw:SmithWaterman = new SmithWaterman(fasta1, fasta2, match, openPanalty, extPanalty);

        x10.Console.OUT.println("The max alignment score: ");
        x10.Console.OUT.println(sw.getMaxScore());

        x10.Console.OUT.println("Matches: ");
        x10.Console.OUT.println(getMatchNumber());

    }

}







class FastaReader {
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






