/*
 * Copyright 2016-2018 Lauri Himanen, Fawzi Mohamed
 * 
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 */

package eu.nomad_lab.parsers

import eu.{ nomad_lab => lab }
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import scala.collection.breakOut

object NWChemParser extends SimpleExternalParserGenerator(
  name = "NWChemParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("NWChemParser")) ::
      ("parserId" -> jn.JString("NWChemParser" + lab.NwchemVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.NwchemVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """              Northwest Computational Chemistry Package \(NWChem\) \d+\.\d+
              ------------------------------------------------------


                    Environmental Molecular Sciences Laboratory
                       Pacific Northwest National Laboratory
                                Richland, WA 99352""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/nwchem/parser/parser-nwchem/nwchemparser/scalainterface.py",
    "${mainFilePath}"),
  cmdCwd = "${mainFilePath}/..",
  resList = Seq(
    "parser-nwchem/nwchemparser/__init__.py",
    "parser-nwchem/nwchemparser/setup_paths.py",
    "parser-nwchem/nwchemparser/parser.py",
    "parser-nwchem/nwchemparser/scalainterface.py",
    "parser-nwchem/nwchemparser/versions/__init__.py",
    "parser-nwchem/nwchemparser/versions/nwchem66/__init__.py",
    "parser-nwchem/nwchemparser/versions/nwchem66/mainparser.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta.nomadmetainfo.json",
    "nomad_meta_info/nwchem.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-nwchem" -> "parsers/nwchem/parser/parser-nwchem",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ DefaultPythonInterpreter.commonDirMapping()
)
