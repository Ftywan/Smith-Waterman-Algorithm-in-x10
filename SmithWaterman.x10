

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
    private var sizesOfVerticalGaps:Array[Short];
    private var sizesOfHorizontalGaps:Array[Short];
    private var prevCells:Array_2[Double];

    public def SmithWatermans(seq1:Rail[String], seq2:Rail[String]) {
        this.seq1 = seq1;
        this.seq2 = seq2;
        this.length1 = seq1.length();
        this.length2 = seq2.length();
        score = new Array_2[Double](length1 + 1, length2 + 1);
        prevCells = new Array_2[Int](length1 + 1, length2 + 1);
   
    }

    publid def similarity(i:Int, j:Int) {
        if (i == 0 || j == 0) {

        }
    }

    public def initMatrix() {
        var i:Int = 0;
        var j:Int = 0;
        for (i = 0; i <= length1; i++) {
            score[i][0] = 0;
            prevCells[i][0] = DR_ZERO;
        }
        for (j = 0; j <= length2; j++) {
            score[0][j] = 0;
            prevCells[0][j] = DR_ZERO;
        }

    }
    public def readFastaFile(fastaFileName:String) {
        val fastaFile = new File(fastaFileName);
        val fastaReader = new FileReader(fastaFile);
        val header = fastaReader.readLine();
        val builder = new StringBuilder();
        val line = fastaReader.readLine();
        while (line != null) {
            builder.append(line);
        }
        return builder.toString();
    }
}