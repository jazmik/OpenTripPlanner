package org.opentripplanner.apis.gtfs.datafetchers;

import graphql.schema.DataFetcher;
import graphql.schema.DataFetchingEnvironment;
import org.locationtech.jts.geom.LineString;
import org.opentripplanner.apis.gtfs.GraphQLRequestContext;
import org.opentripplanner.apis.gtfs.generated.GraphQLDataFetchers;
import org.opentripplanner.apis.gtfs.model.VehiclePositionGeometryModel;
import org.opentripplanner.framework.geometry.SphericalDistanceLibrary;
import org.opentripplanner.framework.geometry.SplitLineString;
import org.opentripplanner.service.realtimevehicles.model.RealtimeVehicle;
import org.opentripplanner.service.realtimevehicles.model.RealtimeVehicle.StopRelationship;
import org.opentripplanner.transit.model.network.TripPattern;
import org.opentripplanner.transit.model.timetable.Trip;
import org.opentripplanner.transit.service.TransitService;

public class VehiclePositionImpl implements GraphQLDataFetchers.GraphQLVehiclePosition {

  @Override
  public DataFetcher<Double> heading() {
    return env -> getSource(env).heading().orElse(null);
  }

  @Override
  public DataFetcher<String> label() {
    return env -> getSource(env).label().orElse(null);
  }

  @Override
  public DataFetcher<Long> lastUpdated() {
    return env -> getSource(env).time().map(time -> time.getEpochSecond()).orElse(null);
  }

  @Override
  public DataFetcher<Double> lat() {
    return env ->
      getSource(env).coordinates().map(coordinates -> coordinates.latitude()).orElse(null);
  }

  @Override
  public DataFetcher<Double> lon() {
    return env ->
      getSource(env).coordinates().map(coordinates -> coordinates.longitude()).orElse(null);
  }

  @Override
  public DataFetcher<Double> odometer()  {
    return env -> getSource(env).odometer().orElse(null);
  }

  @Override
  public DataFetcher<Double> speed() {
    return env -> getSource(env).speed().orElse(null);
  }

  @Override
  public DataFetcher<StopRelationship> stopRelationship() {
    return env -> getSource(env).stop().orElse(null);
  }

  @Override
  public DataFetcher<Trip> trip() {
    return env -> getSource(env).trip();
  }

  @Override
  public DataFetcher<String> vehicleId() {
    return env -> getSource(env).vehicleId().map(vehicleId -> vehicleId.toString()).orElse(null);
  }

  @Override
  public DataFetcher<Object> vehiclePositionGeometry() {
    return environment -> {
      TripPattern pattern = getTripPattern(environment);

      if (pattern == null) {
        return null;
      }

      LineString geometry = pattern.getGeometry();

      if (geometry == null) {
        return null;
      }

      double patternLength = SphericalDistanceLibrary.length(geometry);

      if (patternLength == 0) {
        return null;
      }

      SplitLineString splitPattern = SphericalDistanceLibrary.splitLineStringAtDistance(geometry, getSource(environment).odometer().orElse(0.0));

      return new VehiclePositionGeometryModel(
        splitPattern.beginning().getCoordinates(),
        splitPattern.ending().getCoordinates()
      );
    };
  }

  private RealtimeVehicle getSource(DataFetchingEnvironment environment) {
    return environment.getSource();
  }

  private TripPattern getTripPattern(DataFetchingEnvironment environment) {

    return getTransitService(environment).getPatternForTrip(getSource(environment).trip());
  }

  private TransitService getTransitService(DataFetchingEnvironment environment) {
    return environment.<GraphQLRequestContext>getContext().transitService();
  }



}


