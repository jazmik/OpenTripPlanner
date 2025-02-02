package org.opentripplanner.apis.transmodel.support;

import graphql.ErrorClassification;
import graphql.ExecutionResult;
import graphql.GraphQLError;
import jakarta.ws.rs.core.Response;
import org.opentripplanner.api.json.GraphQLResponseSerializer;
import org.opentripplanner.framework.application.OTPRequestTimeoutException;
import org.opentripplanner.framework.http.OtpHttpStatus;

/**
 * Map the GraphQl execution result to a jakarta Response.
 */
public class ExecutionResultMapper {

  private static final ErrorClassification API_PROCESSING_TIMEOUT = ErrorClassification.errorClassification(
    "ApiProcessingTimeout"
  );
  private static final ErrorClassification BAD_REQUEST_ERROR = ErrorClassification.errorClassification(
    "BadRequestError"
  );

  private static final ErrorClassification INTERNAL_SERVER_ERROR = ErrorClassification.errorClassification(
    "InternalServerError"
  );

  public static Response okResponse(ExecutionResult result) {
    return Response.ok(GraphQLResponseSerializer.serialize(result)).build();
  }

  public static Response timeoutResponse() {
    var error = GraphQLError
      .newError()
      .errorType(API_PROCESSING_TIMEOUT)
      .message(OTPRequestTimeoutException.MESSAGE)
      .build();
    var result = ExecutionResult.newExecutionResult().addError(error).build();
    return response(result, OtpHttpStatus.STATUS_UNPROCESSABLE_ENTITY);
  }

  public static Response badRequestResponse(String message) {
    var error = GraphQLError.newError().errorType(BAD_REQUEST_ERROR).message(message).build();
    var result = ExecutionResult.newExecutionResult().addError(error).build();
    return response(result, Response.Status.BAD_REQUEST);
  }

  public static Response systemErrorResponse(String message) {
    var error = GraphQLError.newError().errorType(INTERNAL_SERVER_ERROR).message(message).build();
    var result = ExecutionResult.newExecutionResult().addError(error).build();
    return response(result, Response.Status.INTERNAL_SERVER_ERROR);
  }

  public static Response response(ExecutionResult result, Response.StatusType status) {
    return Response
      .status(status.getStatusCode())
      .entity(GraphQLResponseSerializer.serialize(result))
      .build();
  }
}
