/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.extensions.picard

import java.io.File
import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction
import htsjdk.samtools.ValidationStringency
import htsjdk.samtools.SAMFileHeader.SortOrder

/**
 * Wraps a Picard function that operates on BAM files.
 * See http://picard.sourceforge.net/ for more info.
 *
 * Since the various BAM utilities take slightly different arguments
 * some values are optional.
 */
trait PicardBamFunction extends JavaCommandLineFunction {
  var validationStringency = ValidationStringency.SILENT
  var sortOrder = SortOrder.coordinate
  var compressionLevel: Option[Int] = None
  var createIndex: Option[Boolean] = None
  var createMD5: Option[Boolean] = None
  var maxRecordsInRam: Option[Int] = None
  var assumeSorted: Option[Boolean] = None

  protected def inputBams: Seq[File]
  protected def outputBam: File

  abstract override def commandLine = super.commandLine +
    repeat("INPUT=", inputBams, spaceSeparated=false) +
    required("TMP_DIR=" + jobLocalDir) +
    optional("OUTPUT=", outputBam, spaceSeparated=false) +
    optional("COMPRESSION_LEVEL=", compressionLevel, spaceSeparated=false) +
    optional("VALIDATION_STRINGENCY=", validationStringency, spaceSeparated=false) +
    optional("SO=", sortOrder, spaceSeparated=false) +
    optional("MAX_RECORDS_IN_RAM=", maxRecordsInRam, spaceSeparated=false) +
    optional("ASSUME_SORTED=", assumeSorted, spaceSeparated=false) +
    optional("CREATE_INDEX=", createIndex, spaceSeparated=false) +
    optional("CREATE_MD5_FILE=", createMD5, spaceSeparated=false)
}
