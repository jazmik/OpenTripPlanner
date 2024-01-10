package org.opentripplanner.apis.gtfs.datafetchers;

import graphql.schema.DataFetcher;
import graphql.schema.DataFetchingEnvironment;
import java.util.List;
import org.locationtech.jts.geom.Coordinate;
import org.opentripplanner.apis.gtfs.generated.GraphQLDataFetchers;
import org.opentripplanner.apis.gtfs.model.PatternHopGeometryModel;

public class PatternHopGeometry implements GraphQLDataFetchers.GraphQLPatternHopGeometry {


  @Override
  public DataFetcher<Double> distance() {
    return environment -> getSource(environment).getDistance();
  }

  @Override
  public DataFetcher<Iterable<Coordinate>> geometry() {
    return environment -> List.of(getSource(environment).getCoordinates());
  }

  private PatternHopGeometryModel getSource(DataFetchingEnvironment environment) {
    return environment.getSource();
  }
}
