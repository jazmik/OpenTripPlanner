package org.opentripplanner.apis.gtfs.model;

import org.locationtech.jts.geom.Coordinate;

/**
 * Class that contains the pattern geometry for a trip split by the vehicle progress
 */
public class VehiclePositionGeometryModel {



  private final Coordinate[] past;

  private final Coordinate[] future;


  public VehiclePositionGeometryModel(Coordinate[] past,Coordinate[] future) {
    this.past = past;
    this.future = future;
  }

  public Coordinate[] getPast() {
    return past;
  }

  public Coordinate[]  getFuture() {
    return future;
  }
}
