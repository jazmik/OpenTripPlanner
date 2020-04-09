package org.opentripplanner.ext.legacygraphqlapi;

import graphql.schema.DataFetcher;
import graphql.schema.idl.TypeRuntimeWiring;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.Arrays;
import java.util.function.Predicate;
import java.util.stream.Collectors;

class IntrospectionTypeWiring {

  static TypeRuntimeWiring build(Class clazz) throws Exception {
    Object instance = clazz.getConstructor().newInstance();

    return TypeRuntimeWiring
        .newTypeWiring(clazz.getSimpleName()
            .replaceFirst("LegacyGraphQL", "")
            .replaceAll("Impl$", ""))
        .dataFetchers(Arrays
            .stream(clazz.getDeclaredMethods())
            .filter(isMethodPublic)
            .filter(isMethodReturnTypeDataFetcher)
            .collect(Collectors.toMap(Method::getName, method -> {
              try { return (DataFetcher) method.invoke(instance); }
              catch (IllegalAccessException | InvocationTargetException ignored) {}
              return null;
            })))
        .build();
  }

  private static Predicate<Method> isMethodPublic = method -> Modifier.isPublic(method.getModifiers());

  private static Predicate<Method> isMethodReturnTypeDataFetcher = (
      (Predicate<Method>) method -> method.getReturnType().equals(DataFetcher.class)
  ).or(method -> Arrays.asList(method.getReturnType().getInterfaces()).contains(DataFetcher.class));
}