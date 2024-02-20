package org.opentripplanner.service.realtimevehicles;

import java.util.List;
import javax.annotation.Nonnull;
import org.locationtech.jts.geom.Envelope;
import org.opentripplanner.service.realtimevehicles.model.RealtimeVehicle;
import org.opentripplanner.transit.model.network.TripPattern;
import org.opentripplanner.transit.model.timetable.OccupancyStatus;
import org.opentripplanner.transit.model.timetable.Trip;

public interface RealtimeVehicleService {
  /**
   * Get the realtime vehicles for a certain trip pattern. Service contains all the vehicles that
   * exist in input feeds but doesn't store any historical data.
   */
  List<RealtimeVehicle> getRealtimeVehicles(@Nonnull TripPattern pattern);

  /**
   * Get the realtime vehicles within the provided Geo Envelope.
   */
  List<RealtimeVehicle> findRealtimeVehiclesInEnvelope(@Nonnull Envelope envelope);

  /**
   * Get the current RealtimeVehicle for a certain trip.
   */
  RealtimeVehicle getRealtimeVehicle(@Nonnull Trip trip);

  /**
   * Get the latest occupancy status for a certain trip.
   */
  OccupancyStatus getVehicleOccupancyStatus(@Nonnull Trip trip);
}
