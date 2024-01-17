package org.opentripplanner.apis.gtfs.datafetchers;

import graphql.schema.DataFetcher;
import graphql.schema.DataFetchingEnvironment;
import java.util.List;
import org.locationtech.jts.geom.Coordinate;
import org.opentripplanner.apis.gtfs.generated.GraphQLDataFetchers;
import org.opentripplanner.apis.gtfs.model.PatternHopGeometryModel;
import org.opentripplanner.apis.gtfs.model.VehiclePositionGeometryModel;

public class VehiclePositionGeometry implements GraphQLDataFetchers.GraphQLVehiclePositionGeometry {


  @Override
  public DataFetcher<Iterable<Coordinate>> future() {
    return environment -> List.of(getSource(environment).getFuture());
  }

  @Override
  public DataFetcher<Iterable<Coordinate>> past() {
    return environment -> List.of(getSource(environment).getPast());
  }

  private VehiclePositionGeometryModel getSource(DataFetchingEnvironment environment) {
    return environment.getSource();
  }
}
