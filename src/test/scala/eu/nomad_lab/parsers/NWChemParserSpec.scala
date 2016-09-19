package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object NWChemParserSpec extends Specification {
  "NWChemParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(NWChemParser, "parsers/nwchem/test/examples/single_point/output.out", "json-events") must_== ParseResult.ParseSuccess
    }
  }

  "test single_point with json" >> {
    ParserRun.parse(NWChemParser, "parsers/nwchem/test/examples/single_point/output.out", "json") must_== ParseResult.ParseSuccess
  }

  "test geo_opt with json" >> {
    ParserRun.parse(NWChemParser, "parsers/nwchem/test/examples/geo_opt/output.out", "json") must_== ParseResult.ParseSuccess
  }

  "test md with json" >> {
    ParserRun.parse(NWChemParser, "parsers/nwchem/test/examples/md/output.out", "json") must_== ParseResult.ParseSuccess
  }
}
