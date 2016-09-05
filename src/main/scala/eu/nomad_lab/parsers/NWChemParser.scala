package eu.nomad_lab.parsers

import eu.{ nomad_lab => lab }
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import scala.collection.breakOut

object NWChemParser extends SimpleExternalParserGenerator(
  name = "NWChemParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("NWChemParser")) ::
      ("parserId" -> jn.JString("NWChemParser" + lab.NWChemVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.NWChemVersionInfo.toMap.map {
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
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/nwchem.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-nwchem" -> "parsers/nwchem/parser/parser-nwchem",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ DefaultPythonInterpreter.commonDirMapping()
)
